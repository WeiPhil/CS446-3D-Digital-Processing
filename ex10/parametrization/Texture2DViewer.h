//
// Created by ziqwang on 16.11.18.
//

#ifndef EX7_SMOOTHING_TEXTURE2DVIEWER_H
#define EX7_SMOOTHING_TEXTURE2DVIEWER_H
#include <nanogui/opengl.h>
#include <nanogui/glutil.h>
#include <nanogui/screen.h>
#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/popupbutton.h>
#include <nanogui/label.h>
#include <nanogui/button.h>
#include <nanogui/textbox.h>
#include <nanogui/tabwidget.h>
#include "mesh_processing.h"

class Texture2DViewer: public nanogui::Widget
{
public:

    nanogui::GLShader shader2D_;

    mesh_processing::MeshProcessing *mesh_;

public:
    Texture2DViewer(Widget* parent, mesh_processing::MeshProcessing *mesh);

    ~Texture2DViewer();

public:

    void init_shader();

    void refresh_mesh();

public:
    void draw(NVGcontext *ctx) override ;

    void drawBorder(NVGcontext* ctx);
};

#endif //EX7_SMOOTHING_TEXTURE2DVIEWER_H
