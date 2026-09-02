// Dear ImGui: standalone example application for GLFW + Metal, using programmable pipeline
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)

// Learn about Dear ImGui:
// - FAQ                  https://dearimgui.com/faq
// - Getting Started      https://dearimgui.com/getting-started
// - Documentation        https://dearimgui.com/docs (same as your local docs/ folder).
// - Introduction, links and more at the top of imgui.cpp

#include <vector>
#include <algorithm>
#include <cstddef>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_metal.h"
#include "mat.h"
#include <stdio.h>

#define GLFW_INCLUDE_NONE
#define GLFW_EXPOSE_NATIVE_COCOA
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>

#import <Metal/Metal.h>
#import <QuartzCore/QuartzCore.h>

#include "scroll_canvas.h"

#include "sym_assert.h"
#include "arena.h"
#include "linearizer.h"
#include "solver.h"

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

typedef struct {
  i32* camera_indices;
  i32* point_indices;
  f64* pixels;
  i32 num_cameras;
  i32 num_points;
  i32 num_observations;
} bal_problem;

typedef struct {
    i32 col;
    i32 row;
} sparse_cell;

typedef struct {
    float canvas_origin[2];
    float cell_size;
    float padding;
    float display_pos[2];
    float display_size[2];
} sparse_uniforms;

static constexpr NSUInteger kVerticesPerCell = 6;

static_assert(offsetof(sparse_uniforms, display_pos) == 4 * sizeof(float));
static_assert(sizeof(sparse_uniforms) == 8 * sizeof(float));

typedef struct {
    id<MTLRenderCommandEncoder> encoder;
    id<MTLRenderPipelineState> pipeline;
    id<MTLBuffer> cells;
    sparse_uniforms uniforms;
    i32 num_cells;
} sparse_draw;

static id<MTLRenderPipelineState> new_sparse_pipeline(
    id<MTLDevice> device,
    MTLPixelFormat pixel_format)
{
    NSString* source = @
        "#include <metal_stdlib>\n"
        "using namespace metal;\n"
        "struct Cell { int col; int row; };\n"
        "struct Uniforms {\n"
        "  float2 canvas_origin;\n"
        "  float cell_size;\n"
        "  float2 display_pos;\n"
        "  float2 display_size;\n"
        "};\n"
        "vertex float4 sparse_vertex(\n"
        "    uint vertex_id [[vertex_id]],\n"
        "    uint instance_id [[instance_id]],\n"
        "    device const Cell* cells [[buffer(0)]],\n"
        "    constant Uniforms& uniforms [[buffer(1)]]) {\n"
        "  constexpr float2 corners[] = {\n"
        "    {0, 0}, {1, 0}, {1, 1}, {0, 0}, {1, 1}, {0, 1}\n"
        "  };\n"
        "  float2 cell = float2(cells[instance_id].col, cells[instance_id].row);\n"
        "  float2 screen = uniforms.canvas_origin +\n"
        "      (cell + corners[vertex_id]) * uniforms.cell_size;\n"
        "  float2 unit = (screen - uniforms.display_pos) / uniforms.display_size;\n"
        "  return float4(unit.x * 2.0 - 1.0, 1.0 - unit.y * 2.0, 0, 1);\n"
        "}\n"
        "fragment half4 sparse_fragment() {\n"
        "  return half4(1);\n"
        "}\n";

    NSError* error = nil;
    id<MTLLibrary> library = [device newLibraryWithSource:source options:nil error:&error];
    if (library == nil) {
        NSLog(@"Sparse shader compilation failed: %@", error);
        return nil;
    }

    MTLRenderPipelineDescriptor* descriptor = [MTLRenderPipelineDescriptor new];
    descriptor.vertexFunction = [library newFunctionWithName:@"sparse_vertex"];
    descriptor.fragmentFunction = [library newFunctionWithName:@"sparse_fragment"];
    descriptor.colorAttachments[0].pixelFormat = pixel_format;

    id<MTLRenderPipelineState> pipeline =
        [device newRenderPipelineStateWithDescriptor:descriptor error:&error];
    if (pipeline == nil) {
        NSLog(@"Sparse pipeline creation failed: %@", error);
    }

    return pipeline;
}

static void draw_sparse_shader(const ImDrawList*, const ImDrawCmd* command)
{
    const sparse_draw* draw = (const sparse_draw*) command->UserCallbackData;
    ImDrawData* data = ImGui::GetDrawData();
    ImVec2 scale = data->FramebufferScale;
    ImVec2 clip_min(
        (command->ClipRect.x - data->DisplayPos.x) * scale.x,
        (command->ClipRect.y - data->DisplayPos.y) * scale.y
    );
    ImVec2 clip_max(
        (command->ClipRect.z - data->DisplayPos.x) * scale.x,
        (command->ClipRect.w - data->DisplayPos.y) * scale.y
    );
    i32 width = (i32) (data->DisplaySize.x * scale.x);
    i32 height = (i32) (data->DisplaySize.y * scale.y);

    clip_min.x = std::max(clip_min.x, 0.0f);
    clip_min.y = std::max(clip_min.y, 0.0f);
    clip_max.x = std::min(clip_max.x, (float) width);
    clip_max.y = std::min(clip_max.y, (float) height);
    if (clip_max.x <= clip_min.x || clip_max.y <= clip_min.y) {
        return;
    }

    MTLScissorRect scissor = {
        .x = (NSUInteger) clip_min.x,
        .y = (NSUInteger) clip_min.y,
        .width = (NSUInteger) (clip_max.x - clip_min.x),
        .height = (NSUInteger) (clip_max.y - clip_min.y),
    };

    [draw->encoder setScissorRect:scissor];
    [draw->encoder setRenderPipelineState:draw->pipeline];
    [draw->encoder setVertexBuffer:draw->cells offset:0 atIndex:0];
    [draw->encoder setVertexBytes:&draw->uniforms
                           length:sizeof(draw->uniforms)
                          atIndex:1];
    [draw->encoder drawPrimitives:MTLPrimitiveTypeTriangle
                      vertexStart:0
                      vertexCount:kVerticesPerCell
                    instanceCount:(NSUInteger) draw->num_cells];
}

int main(int argc, char** argv)
{
    SYM_ASSERT(argc == 2);

    size n = 1L << 30; // 1 GiB seems fine
    printf("Allocating %ld bytes\n", n);
    u8* buf = (u8*) malloc(n);

    sym_arena arena = {
        .beg = buf,
        .end = buf + n,
    };

    sym_allocator alloc_struct = {
        .malloc = sym_arena_malloc,
        .free = sym_arena_free,
        .ctx = &arena
    };
    sym_allocator* alloc = &alloc_struct;

    FILE* file = fopen(argv[1], "r");

    bal_problem p = {};
    fscanf(file, "%d", &p.num_cameras);
    fscanf(file, "%d", &p.num_points);
    fscanf(file, "%d", &p.num_observations);

    p.camera_indices = (i32*) alloc->malloc(p.num_observations * sizeof(i32), alloc->ctx);
    p.point_indices = (i32*) alloc->malloc(p.num_observations * sizeof(i32), alloc->ctx);

    for (i32 i = 0; i < p.num_observations; i++) {
        i32 camera, point;
        fscanf(file, "%d", &camera);
        fscanf(file, "%d", &point);

        f64 px, py;
        fscanf(file, "%lf", &px);
        fscanf(file, "%lf", &py);

        p.camera_indices[i] = camera;
        p.point_indices[i] = point;
    }

    fclose(file);

    // Compute Hessian_lower block triplets.
    i32 nblocks = p.num_observations * 6;
    i32* Hl_block_rows = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
    i32* Hl_block_cols = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);

    // Linearizer callers are responsible for choosing the order of the keys.
    // Here we have all the cameras in order (pose then intrinsics) followed by all the points.
    i32 nkeys = 2 * p.num_cameras + p.num_points;
    for (i32 obs_index = 0; obs_index < p.num_observations; ++obs_index) {
        i32 camera_index = p.camera_indices[obs_index];
        i32 pose_key = 2 * camera_index + 0;
        i32 intrinsics_key = 2 * camera_index + 1;
        i32 point_index = p.point_indices[obs_index];
        i32 point_key = 2 * p.num_cameras + point_index;

        // NOTE: The order here must be consistent with the order later on in the update calls.
        Hl_block_rows[6 * obs_index + 0] = pose_key;
        Hl_block_rows[6 * obs_index + 1] = intrinsics_key;
        Hl_block_rows[6 * obs_index + 2] = point_key;
        Hl_block_rows[6 * obs_index + 3] = intrinsics_key;
        Hl_block_rows[6 * obs_index + 4] = point_key;
        Hl_block_rows[6 * obs_index + 5] = point_key;

        Hl_block_cols[6 * obs_index + 0] = pose_key;
        Hl_block_cols[6 * obs_index + 1] = pose_key;
        Hl_block_cols[6 * obs_index + 2] = pose_key;
        Hl_block_cols[6 * obs_index + 3] = intrinsics_key;
        Hl_block_cols[6 * obs_index + 4] = intrinsics_key;
        Hl_block_cols[6 * obs_index + 5] = point_key;
    }

    // Compute key sizes.
    i32* key_sizes = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
    for (i32 i = 0; i < p.num_cameras; ++i) {
        key_sizes[2 * i + 0] = 6;
        key_sizes[2 * i + 1] = 3;
    }
    for (i32 i = 0; i < p.num_points; ++i) {
        key_sizes[2 * p.num_cameras + i] = 3;
    }

    std::vector<f32> key_sizes_f32;
    for (i32 i = 0; i < nkeys; ++i) {
        key_sizes_f32.push_back(key_sizes[i]);
    }

    // Create the linearizer, linearization.
    i32* Hl_block_nz_indices = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
    sym_csc_mat Hl_block = sym_csc_from_pairs(Hl_block_rows, Hl_block_cols, nblocks, nkeys, nkeys, Hl_block_nz_indices, alloc);
    alloc->free(Hl_block_rows, nblocks * sizeof(i32), alloc->ctx);
    alloc->free(Hl_block_cols, nblocks * sizeof(i32), alloc->ctx);

    i32* key_perm = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
    sym_get_metis_tri_perm(Hl_block, key_sizes, NULL, key_perm, alloc);

    sym_linearization lin;
    sym_linearizer lzr = sym_linearizer_new(
        Hl_block, Hl_block_nz_indices, nblocks,
        key_sizes, nkeys,
        key_perm,
        &lin,
        alloc
    );
    lin.Hl.data = (f64*) alloc->malloc(lin.Hl.nnz * sizeof(f64), alloc->ctx);
    lin.rhs.data = (f64*) alloc->malloc(lin.rhs.n * sizeof(f64), alloc->ctx);

    alloc->free(key_perm, nkeys * sizeof(i32), alloc->ctx);
    alloc->free(key_sizes, nkeys * sizeof(i32), alloc->ctx);

    i32* Hlt_perm = (i32*) alloc->malloc(lin.Hl.nnz * sizeof(i32), alloc->ctx);
    sym_csc_mat Hlt = sym_transpose_csc(lin.Hl, Hlt_perm, alloc);

    sym_chol_factorization fac = {};
    sym_chol_solver solver;
    solver = sym_new_chol_solver(Hlt, &fac, false, alloc);

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    // Create window with graphics context
    float main_scale = ImGui_ImplGlfw_GetContentScaleForMonitor(glfwGetPrimaryMonitor()); // Valid on GLFW 3.3+ only
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    GLFWwindow* window = glfwCreateWindow((int)(1280 * main_scale), (int)(800 * main_scale), "Dear ImGui GLFW+Metal example", nullptr, nullptr);
    if (window == nullptr)
        return 1;

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // Setup scaling
    ImGuiStyle& style = ImGui::GetStyle();
    style.ScaleAllSizes(main_scale);        // Bake a fixed style scale. (until we have a solution for dynamic style scaling, changing this requires resetting Style + calling this again)
    style.FontScaleDpi = main_scale;        // Set initial font scale. (in docking branch: using io.ConfigDpiScaleFonts=true automatically overrides this for every window depending on the current monitor)

    id <MTLDevice> device = MTLCreateSystemDefaultDevice();
    id <MTLCommandQueue> commandQueue = [device newCommandQueue];

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplMetal_Init(device);

    // Load Fonts
    // - If fonts are not explicitly loaded, Dear ImGui will select an embedded font: either AddFontDefaultVector() or AddFontDefaultBitmap().
    //   This selection is based on (style.FontSizeBase * style.FontScaleMain * style.FontScaleDpi) reaching a small threshold.
    // - You can load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - If a file cannot be loaded, AddFont functions will return a nullptr. Please handle those errors in your code (e.g. use an assertion, display an error and quit).
    // - Read 'docs/FONTS.md' for more instructions and details.
    // - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use FreeType for higher quality font rendering.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //style.FontSizeBase = 20.0f;
    //io.Fonts->AddFontDefaultVector();
    //io.Fonts->AddFontDefaultBitmap();
    //io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf");
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf");
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf");
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf");
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf");
    ImFont* font = io.Fonts->AddFontFromFileTTF("/System/Library/Fonts/SFNSMono.ttf", 16.0f);
    IM_ASSERT(font != nullptr);

    NSWindow *nswin = glfwGetCocoaWindow(window);
    CAMetalLayer *layer = [CAMetalLayer layer];
    layer.device = device;
    layer.pixelFormat = MTLPixelFormatBGRA8Unorm;
    nswin.contentView.layer = layer;
    nswin.contentView.wantsLayer = YES;

    MTLRenderPassDescriptor *renderPassDescriptor = [MTLRenderPassDescriptor new];

    id<MTLRenderPipelineState> sparse_pipeline =
        new_sparse_pipeline(device, layer.pixelFormat);
    SYM_ASSERT(sparse_pipeline != nil);

    // Our state
    float clear_color[4] = {0.45f, 0.55f, 0.60f, 1.00f};

    // Main loop
    while (!glfwWindowShouldClose(window))
    {
        @autoreleasepool
        {
            // Poll and handle events (inputs, window resize, etc.)
            // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
            // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
            // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
            // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
            glfwPollEvents();

            int width, height;
            glfwGetFramebufferSize(window, &width, &height);
            layer.drawableSize = CGSizeMake(width, height);
            id<CAMetalDrawable> drawable = [layer nextDrawable];

            id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
            renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(clear_color[0] * clear_color[3], clear_color[1] * clear_color[3], clear_color[2] * clear_color[3], clear_color[3]);
            renderPassDescriptor.colorAttachments[0].texture = drawable.texture;
            renderPassDescriptor.colorAttachments[0].loadAction = MTLLoadActionClear;
            renderPassDescriptor.colorAttachments[0].storeAction = MTLStoreActionStore;
            id <MTLRenderCommandEncoder> renderEncoder = [commandBuffer renderCommandEncoderWithDescriptor:renderPassDescriptor];
            [renderEncoder pushDebugGroup:@"ImGui demo"];

            // Start the Dear ImGui frame
            ImGui_ImplMetal_NewFrame(renderPassDescriptor);
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            const ImGuiViewport *vp = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(vp->WorkPos);
            ImGui::SetNextWindowSize(vp->WorkSize);

            /* Pinned to the viewport, so the decorations a floating window would need
            are only in the way. */
            const ImGuiWindowFlags flags =
                ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus |
                ImGuiWindowFlags_NoSavedSettings;

            i32 col_starts[] = {0, 1, 4, 5, 6, 7};
            i32 row_indices[] = {0, 1, 3, 4, 2, 3, 4};

            sym_csc_mat m = {
                .col_starts = col_starts,
                .row_indices = row_indices,
                .data = nullptr,
                .nrows = 5,
                .ncols = 5,
                .nnz = 7,
            };

            if (ImGui::Begin("ImGui Base", nullptr, flags)) {
                static bool use_sparse_shader = false;
                ImGui::Checkbox("Dedicated shader", &use_sparse_shader);

                // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
                ImVec2 canvas_min = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
                ImVec2 canvas_size = ImGui::GetContentRegionAvail();
                ImVec2 canvas_max = ImVec2(canvas_min.x + canvas_size.x, canvas_min.y + canvas_size.y);

                ImGui::InvisibleButton("canvas", canvas_size, ImGuiButtonFlags_MouseButtonLeft);

                ImDrawList* draw_list = ImGui::GetWindowDrawList();

                static ScrollCanvas canvas;

                const auto draw_mat = [&draw_list](sym_csc_mat m, f32 base_size, f32* row_weights = nullptr, f32* col_weights = nullptr, bool fill_diag = false) {
                    std::vector<f32> row_size_scan;
                    row_size_scan.push_back(0.0f);
                    for (i32 i = 0; i < m.nrows; ++i) {
                        row_size_scan.push_back(row_size_scan.back() + (row_weights == nullptr ? base_size : row_weights[i] * base_size));
                    }

                    std::vector<f32> col_size_scan;
                    col_size_scan.push_back(0.0f);
                    for (i32 i = 0; i < m.ncols; ++i) {
                        col_size_scan.push_back(col_size_scan.back() + (col_weights == nullptr ? base_size : col_weights[i] * base_size));
                    }

                    draw_list->AddRectFilled(
                        canvas.CanvasToViewport(ImVec2(col_size_scan.front(), row_size_scan.front())),
                        canvas.CanvasToViewport(ImVec2(col_size_scan.back(), row_size_scan.back())),
                        IM_COL32(25, 25, 25, 255)
                    );

                    // for (i32 r = 0; r <= m.nrows; ++r) {
                    //     draw_list->AddLine(
                    //         canvas.CanvasToViewport(ImVec2(0.0, row_size_scan.at(r))),
                    //         canvas.CanvasToViewport(ImVec2(col_size_scan.back(), row_size_scan.at(r))),
                    //         IM_COL32(150, 150, 150, 255)
                    //     );
                    // }
                    // for (i32 c = 0; c <= m.ncols; ++c) {
                    //     draw_list->AddLine(
                    //         canvas.CanvasToViewport(ImVec2(col_size_scan.at(c), 0.0)),
                    //         canvas.CanvasToViewport(ImVec2(col_size_scan.at(c), row_size_scan.back())),
                    //         IM_COL32(150, 150, 150, 255)
                    //     );
                    // }

                    int c = 0;
                    for (i32 i = 0; i < m.nnz; ++i) {
                        while (m.col_starts[c + 1] <= i) {
                            ++c;
                        }
                        int r = m.row_indices[i];
                        draw_list->AddRectFilled(
                            canvas.CanvasToViewport(ImVec2(col_size_scan.at(c), row_size_scan.at(r))),
                            canvas.CanvasToViewport(ImVec2(col_size_scan.at(c + 1), row_size_scan.at(r + 1))),
                            IM_COL32(255, 255, 255, 255)
                        );
                    }

                    if (fill_diag) {
                        for (i32 i = 0; i < m.nrows && i < m.ncols; ++i) {
                            draw_list->AddRectFilled(
                                canvas.CanvasToViewport(ImVec2(col_size_scan.at(i), row_size_scan.at(i))),
                                canvas.CanvasToViewport(ImVec2(col_size_scan.at(i + 1), row_size_scan.at(i + 1))),
                                IM_COL32(255, 255, 255, 255)
                            );
                        }
                    }
                };

                const auto draw_mat_shader = [&](sym_csc_mat m, f32 base_size, bool fill_diag = false) {
                    ImVec2 matrix_min = canvas.CanvasToViewport(ImVec2(0.0f, 0.0f));
                    ImVec2 matrix_max = canvas.CanvasToViewport(ImVec2(
                        m.ncols * base_size,
                        m.nrows * base_size
                    ));
                    draw_list->AddRectFilled(
                        matrix_min,
                        matrix_max,
                        IM_COL32(25, 25, 25, 255)
                    );

                    // One shader instance draws one CSC nonzero.
                    std::vector<sparse_cell> cells;
                    cells.reserve(m.nnz);
                    for (i32 col = 0; col < m.ncols; ++col) {
                        for (i32 nz = m.col_starts[col]; nz < m.col_starts[col + 1]; ++nz) {
                            cells.push_back({col, m.row_indices[nz]});
                        }
                    }
                    if (fill_diag) {
                        for (i32 i = 0; i < m.nrows && i < m.ncols; ++i) {
                            cells.push_back({i, i});
                        }
                    }
                    if (cells.empty()) {
                        return;
                    }

                    id<MTLBuffer> buffer = [[device
                        newBufferWithBytes:cells.data()
                                    length:cells.size() * sizeof(sparse_cell)
                                   options:MTLResourceStorageModeShared] autorelease];
                    sparse_draw draw = {
                        .encoder = renderEncoder,
                        .pipeline = sparse_pipeline,
                        .cells = buffer,
                        .uniforms = {
                            .canvas_origin = {matrix_min.x, matrix_min.y},
                            .cell_size = base_size * canvas.Scale(),
                            .display_pos = {vp->Pos.x, vp->Pos.y},
                            .display_size = {vp->Size.x, vp->Size.y},
                        },
                        .num_cells = (i32) cells.size(),
                    };
                    draw_list->AddCallback(draw_sparse_shader, &draw, sizeof(draw));
                    draw_list->AddCallback(ImDrawCallback_ResetRenderState, nullptr);
                };

                {
                    canvas.Begin(canvas_min, canvas_max, {});

                    draw_list->PushClipRect(canvas_min, canvas_max, true);

                    draw_list->AddRectFilled(canvas_min, canvas_max, IM_COL32(50, 50, 50, 255));

                    // std::vector<f32> row_weights(m.nrows, 1.0);
                    // row_weights.at(1) = 2.0;
                    // std::vector<f32> col_weights(m.ncols, 1.0);
                    // col_weights.at(4) = 2.0;
                    constexpr auto kCellSize = 50.0f;
                    // draw_mat(m, kCellSize, row_weights.data(), col_weights.data());


                    // draw_mat(Hl_block, kCellSize, key_sizes_f32.data(), key_sizes_f32.data());
                    // draw_mat(Hl_block, kCellSize);
                    if (use_sparse_shader) {
                        draw_mat_shader(fac.L, kCellSize, true);
                    } else {
                        draw_mat(fac.L, kCellSize, nullptr, nullptr, true);
                    }

                    canvas.End();
                }

                ImGui::End();
            }

            // Rendering
            ImGui::Render();
            ImGui_ImplMetal_RenderDrawData(ImGui::GetDrawData(), commandBuffer, renderEncoder);

            [renderEncoder popDebugGroup];
            [renderEncoder endEncoding];

            [commandBuffer presentDrawable:drawable];
            [commandBuffer commit];
        }
    }

    // Cleanup
    ImGui_ImplMetal_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
