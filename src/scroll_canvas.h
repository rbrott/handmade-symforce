#pragma once

#include <optional>
#include <algorithm>
#include <cmath>

#include "imgui.h"

inline bool InRect(const ImVec2& rect_min, const ImVec2& rect_max,
                   const ImVec2& query) {
    return rect_min.x <= query.x && rect_min.y <= query.y && query.x <= rect_max.x && query.y <= rect_max.y;
}

class ScrollCanvas {
  public:
    struct Options {
        // Temporarily override the normal user input.
        std::optional<ImVec2> center_override;
        std::optional<double> scale_override;

        bool add_grid{false};
    };

    void Begin(const ImVec2& viewport_min, const ImVec2& viewport_max, Options options);
    void End();

    ImVec2 CanvasToViewport(const ImVec2& p) const;
    ImVec2 ViewportToCanvas(const ImVec2& p) const;

    bool ContainsViewport(const ImVec2& p) const {
        return InRect(viewport_min_, viewport_max_, p);
    }

    float Scale() const {
        return canvas_scale_;
    }

  private:
    ImVec2 viewport_min_{};  // in screen coords
    ImVec2 viewport_max_{};  // in screen coords

    ImVec2 canvas_origin_{};    // in viewport coords
    float canvas_scale_{1.0f};

    ImVec2 canvas_origin_input_{};    // in viewport coords
    float canvas_scale_input_{1.0f};
};

void ScrollCanvas::Begin(const ImVec2& viewport_min, const ImVec2& viewport_max, Options options) {
    ImGuiIO& io = ImGui::GetIO();

    if (options.scale_override.has_value()) {
        canvas_scale_ = options.scale_override.value();
    } else {
        if (ImGui::IsWindowHovered() && io.MouseWheel != 0.0f && !ImGui::IsAnyItemActive() && InRect(viewport_min, viewport_max, io.MousePos)) {
            double new_scale = canvas_scale_input_ * std::pow(1.2, 0.25 * io.MouseWheel);
            // TODO: Clamping here vs. widget?
            // new_scale = std::clamp(new_scale, 0.1, 50.0);

            canvas_origin_input_.x += (canvas_scale_input_ - new_scale) * ViewportToCanvas(io.MousePos).x;
            canvas_origin_input_.y += (canvas_scale_input_ - new_scale) * ViewportToCanvas(io.MousePos).y;
            canvas_scale_input_ = new_scale;
        }

        canvas_origin_ = canvas_origin_input_;
        canvas_scale_ = canvas_scale_input_;
    }

    if (options.center_override.has_value()) {
        canvas_origin_.x = 0.5 * (viewport_max.x - viewport_min.x) - canvas_scale_ * options.center_override.value().x;
        canvas_origin_.y = 0.5 * (viewport_max.y - viewport_min.y) - canvas_scale_ * options.center_override.value().y;
    } else {
        if (ImGui::IsItemActive() && ImGui::IsMouseDragging(ImGuiMouseButton_Left, /*threshold=*/0.0f) && InRect(viewport_min, viewport_max, io.MousePos)) {
            canvas_origin_input_.x += ImVec2(io.MouseDelta).x;
            canvas_origin_input_.y += ImVec2(io.MouseDelta).y;
        }

        canvas_origin_ = canvas_origin_input_;
    }

    viewport_min_ = viewport_min;
    viewport_max_ = viewport_max;
}

void ScrollCanvas::End() {}

ImVec2 ScrollCanvas::CanvasToViewport(const ImVec2& p) const {
    return ImVec2(
        viewport_min_.x + canvas_origin_.x + canvas_scale_ * p.x,
        viewport_min_.y + canvas_origin_.y + canvas_scale_ * p.y
    );
}

ImVec2 ScrollCanvas::ViewportToCanvas(const ImVec2& p) const {
    return ImVec2(
        (p.x - viewport_min_.x - canvas_origin_.x) / canvas_scale_,
        (p.y - viewport_min_.y - canvas_origin_.y) / canvas_scale_
    );
}

