// Dear ImGui: standalone example application for GLFW + Metal, using programmable pipeline
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)

// Learn about Dear ImGui:
// - FAQ                  https://dearimgui.com/faq
// - Getting Started      https://dearimgui.com/getting-started
// - Documentation        https://dearimgui.com/docs (same as your local docs/ folder).
// - Introduction, links and more at the top of imgui.cpp

#include <vector>

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
                // Using InvisibleButton() as a convenience 1) it will advance the layout cursor and 2) allows us to use IsItemHovered()/IsItemActive()
                ImVec2 canvas_min = ImGui::GetCursorScreenPos();      // ImDrawList API uses screen coordinates!
                ImVec2 canvas_size = ImGui::GetContentRegionAvail();
                ImVec2 canvas_max = ImVec2(canvas_min.x + canvas_size.x, canvas_min.y + canvas_size.y);

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
                    draw_mat(fac.L, kCellSize, nullptr, nullptr, true);

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
