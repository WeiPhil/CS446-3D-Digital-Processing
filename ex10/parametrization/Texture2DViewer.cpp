//
// Created by ziqwang on 16.11.18.
//

#include "Texture2DViewer.h"

using Eigen::Vector2i;
using Eigen::Vector2f;
using nanogui::Screen;

Texture2DViewer::Texture2DViewer(nanogui::Widget *parent, mesh_processing::MeshProcessing *mesh)
        :Widget(parent)
{
    mesh_ = mesh;
    mPos = Vector2i(25, 35);
    mSize = Vector2i(350, 350);

    init_shader();
}

Texture2DViewer::~Texture2DViewer() {
    shader2D_.free();
}

void Texture2DViewer::draw(NVGcontext *ctx) {

    Widget::draw(ctx);
    nvgEndFrame(ctx);
    drawBorder(ctx);



    const Screen* screen = dynamic_cast<const Screen*>(this->window()->parent());
    assert(screen);
    float r = screen->pixelRatio();
    Vector2f screenSize = screen->size().cast<float>();
    Vector2f positionInScreen = absolutePosition().cast<float>();
    Vector2f positionAfterOffset = positionInScreen + Vector2f(0, mSize.y()) - screenSize / 2;
    Vector2f origin = positionAfterOffset.cwiseQuotient(screenSize); origin = origin * r;
    Vector2f ratio = Vector2f(mSize.x(), mSize.y()); ratio = ratio.cwiseQuotient(screenSize / 2);
    glEnable(GL_SCISSOR_TEST);
    glScissor(positionInScreen.x() * r,
              (screenSize.y() - positionInScreen.y() - size().y()) * r,
              size().x() * r, size().y() * r);


    // Draw Lines;
    shader2D_.bind();
    shader2D_.setUniform("origin", origin);
    shader2D_.setUniform("ratio", ratio);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    shader2D_.drawIndexed(GL_TRIANGLES, 0, mesh_->get_number_of_face());
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glDisable(GL_SCISSOR_TEST);
}


void Texture2DViewer::refresh_mesh()
{
    shader2D_.bind();
    shader2D_.uploadIndices(*(mesh_->get_indices()));
    shader2D_.uploadAttrib("vertex", *(mesh_->get_textures()));
}


void Texture2DViewer::drawBorder(NVGcontext *ctx)
{
    nvgSave(ctx);
    nvgBeginPath(ctx);
    nvgScissor(ctx, mPos.x(), mPos.y(), mSize.x(), mSize.y());
    nvgStrokeWidth(ctx, 3.0f);
    Vector2i borderPosition = mPos;
    Vector2i borderSize = mSize;
    nvgRect(ctx, borderPosition.x() - 0.5f, borderPosition.y() - 0.5f,
            borderSize.x() + 1, borderSize.y() + 1);
    nvgStrokeColor(ctx, nanogui::Color(0.0f, 0.0f, 0.0f, 1.0f));
    nvgStroke(ctx);
    nvgResetScissor(ctx);
    nvgRestore(ctx);
}

void Texture2DViewer::init_shader()
{
    shader2D_.init(
            /* An identifying name */
            "a_simple_shader",

            /* Vertex shader */
            "#version 330\n"
            "in vec2 vertex;\n"
            "uniform vec2 origin;\n"
            "uniform vec2 ratio;\n"
            "void main() {\n"
            "    gl_Position = vec4(vertex.x * ratio.x + origin.x, vertex.y * ratio.y - origin.y, 0.0, 1.0);\n"
            "}",

            /* Fragment shader */
            "#version 330\n"
            "out vec4 color;\n"
            "void main() {\n"
            " color = vec4(1.0, 1.0, 1.0, 1.0);\n"
            "}"
    );
}
