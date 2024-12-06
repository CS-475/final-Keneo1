/*
 *  Copyright 2024 <me>
 */

#ifndef _g_starter_canvas_h_
#define _g_starter_canvas_h_

#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GBitmap.h"
#include "include/GBlendMode.h"
#include "include/GPaint.h"
#include "include/GMatrix.h"
#include <stack>

class MyCanvas : public GCanvas {
public:
    MyCanvas(const GBitmap& device) : fDevice(device) {
        stackMatrix.push(GMatrix());
    }

    void clear(const GColor& color) override;
    void drawRect(const GRect& rect, const GPaint& color) override;
    void drawConvexPolygon(const GPoint points[], int count, const GPaint& color) override;
    virtual void drawPath(const GPath&, const GPaint&) override;
    virtual void blitrow(GPixel row_addr[], int N, int x, int y, GPixel src, GBlendMode mode, GShader* shader);
    virtual void save() override;
    virtual void restore() override;
    virtual void concat(const GMatrix&) override;
    virtual void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint&) override;
    virtual void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint&) override;

private:
    // Note: we store a copy of the bitmap
    const GBitmap fDevice;
    std::stack<GMatrix> stackMatrix;
    GMatrix CTM = GMatrix();

    // Add whatever other fields you need
};
#endif


