/*
 *  Copyright 2024 <me>
 */

#include "my_canvas.h"
#include "include/GRect.h"
#include "include/GPoint.h"
#include "include/GBlendMode.h"
#include "include/GPath.h"
#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GColor.h"
#include "include/GPathBuilder.h"
#include <assert.h>
#include <iostream>
#include <vector>
#include <stack>
#include "include/GFinal.h"

using namespace std;


float getQPt1(float A, float B, float C, float t) {
  float qpt = A + t * (B - A) + t * ((B + t * (C - B)) - (A + t * (B - A)));

  return qpt;
}

float getCPt1(float A, float B, float C, float D, float t) {
  float quad1 = getQPt1(A, B, C, t);
  float quad2 = getQPt1(B, C, D, t);
  float cubic1 = quad1 + t * (quad2 - quad1);

  return cubic1;
}


struct Edge {
    float m;
    float b;
    int top;
    int bottom;
    int wind;
};

float eval(Edge e, float y) {
  return e.m * (y + 0.5) + e.b;

}

GPixel ToGPixel(GColor color) {

    int a = int(color.a * 255 + 0.5f);
    int r = int(color.r * color.a * 255 + 0.5f);
    int g = int(color.g * color.a * 255 + 0.5f);
    int b = int(color.b * color.a * 255 + 0.5f);
    return GPixel_PackARGB(a, r, g, b);

}


int myRound(float x) {
    if (x < 0 ) {
        return 0;
    }

    return (int) floor(x + 0.5f);
}

void MyCanvas::clear(const GColor& color) {
    // your code here
    GPixel *row_addr = nullptr;
    GPixel pixel_color = ToGPixel(color);
    int width = fDevice.width();
    int height = fDevice.height();
    for (int y=0; y<height; y++){
        row_addr = fDevice.getAddr(0,y);
        for (int x=0; x<width; x++) {
            row_addr[x] = pixel_color;
        }
    }
}

bool compare(Edge first, Edge second) {
    if (first.top == second.top) {
      return eval(first, first.top) < eval(second, second.top);
    }
    return (first.top < second.top);
}

GIRect changeRect(GRect s_rect, int width, int height) {
    GIRect rect;
    float l = s_rect.left;
    float r = s_rect.right;
    float t = s_rect.top;
    float b = s_rect.bottom;
    rect.left = l < 0 ? 0 : (int)(0.5 + l);
    rect.right = r > width ? width : (int)(0.5 + r);
    rect.top = t < 0 ? 0 : (int)(0.5 + t);
    rect.bottom = b > height ? height : (int)(0.5 + b);
    return rect;

}

int divide255(int x) {
    return (x+128)*257 >>16;
}

GPixel blendPix(GPixel dstPixel, GPixel srcPixel) {
    
    int srcA = GPixel_GetA(srcPixel);
    int srcR = GPixel_GetR(srcPixel);
    int srcG = GPixel_GetG(srcPixel);
    int srcB = GPixel_GetB(srcPixel);

    int dstA = GPixel_GetA(dstPixel);
    int dstR = GPixel_GetR(dstPixel);
    int dstG = GPixel_GetG(dstPixel);
    int dstB = GPixel_GetB(dstPixel);

    int outA = srcA + divide255(dstA * (255 - srcA));
    int outR = srcR + divide255(dstR * (255 - srcA));
    int outG = srcG + divide255(dstG * (255 - srcA));
    int outB = srcB + divide255(dstB * (255 - srcA));

    return GPixel_PackARGB(outA, outR, outG, outB);
}

// kClear: Always returns a fully transparent pixel
GPixel kClear_mode(GPixel dst, GPixel src) {
    return GPixel_PackARGB(0, 0, 0, 0);  // Fully transparent
}

// kSrc: Completely replaces the destination pixel with the source
GPixel kSrc_mode(GPixel dst, GPixel src) {
    return src;
}

// kDst: Keeps the destination pixel, ignores the source
GPixel kDst_mode(GPixel dst, GPixel src) {
    return dst;
}
// kSrcOver: Source blends over destination based on source alpha
GPixel kSrcOver_mode(GPixel dst, GPixel src) {
    int srcA = GPixel_GetA(src);
    int dstA = GPixel_GetA(dst);

    int outA = srcA + divide255(dstA * (255 - srcA));
    int outR = GPixel_GetR(src) + divide255(GPixel_GetR(dst) * (255 - srcA));
    int outG = GPixel_GetG(src) + divide255(GPixel_GetG(dst) * (255 - srcA));
    int outB = GPixel_GetB(src) + divide255(GPixel_GetB(dst) * (255 - srcA));

    return GPixel_PackARGB(outA, outR, outG, outB);
}

// kDstOver: Destination blends over source based on destination alpha
GPixel kDstOver_mode(GPixel dst, GPixel src) {
    int dstA = GPixel_GetA(dst);
    int srcA = GPixel_GetA(src);

    int outA = dstA + divide255(srcA * (255 - dstA));
    int outR = GPixel_GetR(dst) + divide255(GPixel_GetR(src) * (255 - dstA));
    int outG = GPixel_GetG(dst) + divide255(GPixel_GetG(src) * (255 - dstA));
    int outB = GPixel_GetB(dst) + divide255(GPixel_GetB(src) * (255 - dstA));

    return GPixel_PackARGB(outA, outR, outG, outB);
}
// kSrcIn: Source pixel is kept where destination alpha is non-zero
GPixel kSrcIn_mode(GPixel dst, GPixel src) {
    int dstA = GPixel_GetA(dst);
    int outA = divide255(dstA * GPixel_GetA(src));
    int outR = divide255(dstA * GPixel_GetR(src));
    int outG = divide255(dstA * GPixel_GetG(src));
    int outB = divide255(dstA * GPixel_GetB(src));
    return GPixel_PackARGB(outA, outR, outG, outB);
}

// kDstIn: Destination pixel is kept where source alpha is non-zero
GPixel kDstIn_mode(GPixel dst, GPixel src) {
    int srcA = GPixel_GetA(src);
    int outA = divide255(srcA * GPixel_GetA(dst));
    int outR = divide255(srcA * GPixel_GetR(dst));
    int outG = divide255(srcA * GPixel_GetG(dst));
    int outB = divide255(srcA * GPixel_GetB(dst));
    return GPixel_PackARGB(outA, outR, outG, outB);
}

// kSrcOut: Source pixel is kept where destination alpha is zero
GPixel kSrcOut_mode(GPixel dst, GPixel src) {
    int dstA = GPixel_GetA(dst);
    int outA = divide255((255 - dstA) * GPixel_GetA(src));
    int outR = divide255((255 - dstA) * GPixel_GetR(src));
    int outG = divide255((255 - dstA) * GPixel_GetG(src));
    int outB = divide255((255 - dstA) * GPixel_GetB(src));
    return GPixel_PackARGB(outA, outR, outG, outB);
}

// kDstOut: Destination pixel is kept where source alpha is zero
GPixel kDstOut_mode(GPixel dst, GPixel src) {
    int srcA = GPixel_GetA(src);
    int outA = divide255((255 - srcA) * GPixel_GetA(dst));
    int outR = divide255((255 - srcA) * GPixel_GetR(dst));
    int outG = divide255((255 - srcA) * GPixel_GetG(dst));
    int outB = divide255((255 - srcA) * GPixel_GetB(dst));
    return GPixel_PackARGB(outA, outR, outG, outB);
}

// kSrcATop: Source is placed atop destination, and destination pixels remain only where source is transparent
GPixel kSrcATop_mode(GPixel dst, GPixel src) {
    int srcA = GPixel_GetA(src);
    int dstA = GPixel_GetA(dst);
    int outA = divide255(dstA * srcA) + divide255((255 - srcA) * dstA);
    int outR = divide255(dstA * GPixel_GetR(src)) + divide255((255 - srcA) * GPixel_GetR(dst));
    int outG = divide255(dstA * GPixel_GetG(src)) + divide255((255 - srcA) * GPixel_GetG(dst));
    int outB = divide255(dstA * GPixel_GetB(src)) + divide255((255 - srcA) * GPixel_GetB(dst));
    return GPixel_PackARGB(outA, outR, outG, outB);
}



// kSrcATop: Source is placed atop destination, and destination pixels remain only where source is transparent
GPixel kDstATop_mode(GPixel dst, GPixel src) {
    int srcA = GPixel_GetA(src);
    int dstA = GPixel_GetA(dst);
    int outA = divide255(srcA * dstA) + divide255((255 - dstA) * srcA);
    int outR = divide255(srcA * GPixel_GetR(dst)) + divide255((255 - dstA) * GPixel_GetR(src));
    int outG = divide255(srcA * GPixel_GetG(dst)) + divide255((255 - dstA) * GPixel_GetG(src));
    int outB = divide255(srcA * GPixel_GetB(dst)) + divide255((255 - dstA) * GPixel_GetB(src));
    return GPixel_PackARGB(outA, outR, outG, outB);
}

// kXor: Source and destination pixels remain only where the other is transparent
GPixel kXor_mode(GPixel dst, GPixel src) {
    int srcA = GPixel_GetA(src);
    int dstA = GPixel_GetA(dst);
    int outA = divide255((255 - srcA) * dstA) + divide255((255 - dstA) * srcA);
    int outR = divide255((255 - srcA) * GPixel_GetR(dst)) + divide255((255 - dstA) * GPixel_GetR(src));
    int outG = divide255((255 - srcA) * GPixel_GetG(dst)) + divide255((255 - dstA) * GPixel_GetG(src));
    int outB = divide255((255 - srcA) * GPixel_GetB(dst)) + divide255((255 - dstA) * GPixel_GetB(src));
    return GPixel_PackARGB(outA, outR, outG, outB);
}


// Main blending function
GPixel blend(GPixel dst, GPixel src, GBlendMode mode) {
    switch (mode) {
        case GBlendMode::kClear:    return kClear_mode(dst, src);
        case GBlendMode::kSrc:      return kSrc_mode(dst, src);
        case GBlendMode::kDst:      return kDst_mode(dst, src);
        case GBlendMode::kSrcOver:  return kSrcOver_mode(dst, src);
        case GBlendMode::kDstOver:  return kDstOver_mode(dst, src);
        case GBlendMode::kSrcIn:    return kSrcIn_mode(dst, src);
        case GBlendMode::kDstIn:    return kDstIn_mode(dst, src);
        case GBlendMode::kSrcOut:   return kSrcOut_mode(dst, src);
        case GBlendMode::kDstOut:   return kDstOut_mode(dst, src);
        case GBlendMode::kSrcATop:  return kSrcATop_mode(dst, src);
        case GBlendMode::kDstATop:  return kDstATop_mode(dst, src);
        case GBlendMode::kXor:      return kXor_mode(dst, src);
        default:                    assert(false);  // Default to src over
    }
    return 0;
}

void MyCanvas::blitrow(GPixel row_addr[], int N, int x, int y, GPixel src, GBlendMode mode, GShader* shader) {
  if (shader == nullptr) {
    for (int i = 0; i < N; i++) {
        GPixel dst = row_addr[i];
        row_addr[i] = blend(dst, src, mode);
    }
  } else {
    if (!shader->setContext(CTM)) {
      return;
    }

    GPixel row[N];
    shader->shadeRow(x, y, N, row);
    for (int i = 0; i < N; i++) {
      GPixel dst = row_addr[i];
      row_addr[i] = blend(dst, row[i], mode);
    }
  }
}

Edge makeEdge(GPoint p0, GPoint p1, int wind) {
    
    if (p0.y > p1.y) {
      swap(p0, p1);
    }
    int w;
    if (wind) {
      w = 1;
    }
    else {
      w = -1;
    }

    Edge e;
    e.top = GRoundToInt(p0.y);
    e.bottom = GRoundToInt(p1.y);

    e.m = (p1.x - p0.x) / (p1.y - p0.y);
    e.b = p0.x - (p0.y * e.m);
    e.wind = w;

    return e;
}

vector<Edge> clipEdges(int bottom, int right, GPoint p0, GPoint p1) {
  vector<Edge> edges;
  bool wind;

  if (p0.y > p1.y) {
    std::swap(p0, p1);
    wind = false;
  } else {
    wind = true;
  }

  if (p0.y == p1.y) {
    return edges;
  }


  if (p1.y <= 0 || p0.y >= bottom) {
    return edges;
  }

  const float m = (p1.x - p0.x) / (p1.y - p0.y);
  const float b = p0.x - (p0.y * m);

  if (p0.y < 0) {
    p0.x = b;
    p0.y = 0;
  }

  if (p1.y > bottom) {
    p1.x = m * bottom + b;
    p1.y = bottom;
  }

  if (p0.x > p1.x) {
    swap(p0, p1);
  }

  if (p1.x <= 0) { 
    Edge edge = makeEdge({0, p0.y}, {0, p1.y}, wind);

    if (edge.top != edge.bottom) {
      edges.push_back(edge);
    }

    return edges;
  }

  if (p0.x >= right) { 
    Edge edge = makeEdge({right, p0.y}, {right, p1.y}, wind);
    
    if (edge.top != edge.bottom) {
      edges.push_back(edge);
    }

    return edges;
  }
  
  if (p0.x < 0) { 
    float vY = (b / m) * -1;

    Edge edge = makeEdge({0, p0.y}, {0, vY}, wind);

    p0 = {0, vY};

    if (edge.top != edge.bottom) {
      edges.push_back(edge);
    }
  }
  
  if (p1.x > right) {
    float vY = (right-b) / m;

    Edge edge = makeEdge({right, p1.y}, {right, vY}, wind);

    p1 = {(float) right, vY};

    if (edge.top != edge.bottom) {
      edges.push_back(edge);
    }
  }

  Edge edge = makeEdge(p0, p1, wind);
    if (edge.top != edge.bottom) {
      edges.push_back(edge);
    }

  return edges;
}

vector<Edge> createEdges(GBitmap fDevice, int count, const GPoint points[]) {

  vector<Edge> edges;

  for (int i=0; i<count; i++) {
    GPoint p0 = points[i];
    GPoint p1 = points[(i+1) % count];

    vector<Edge> new_edges = clipEdges(fDevice.height(), fDevice.width(), p0, p1);

    if (new_edges.size() > 0) {
      edges.insert(edges.end(), new_edges.begin(), new_edges.end());
    }
  }
  return edges;
}


void MyCanvas::drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) {

    GBlendMode mode = paint.getBlendMode();
    GPixel pixel_color = ToGPixel(paint.getColor());
    GShader* shader = paint.peekShader();
    if (mode == GBlendMode::kDst) {
      return;
    }
    int srcA = GPixel_GetA(pixel_color);
    if (!shader && srcA == 255 || shader && shader->isOpaque()) {
        if (mode == GBlendMode::kSrcOver) {
            mode = GBlendMode::kSrc;
        }
        else if (mode == GBlendMode::kDstIn) {
            return;
        }
        else if (mode == GBlendMode::kDstOut) {
            mode = GBlendMode::kClear;
        }
        else if (mode == GBlendMode::kSrcATop) {
            mode = GBlendMode::kSrcIn;
        }
        else if (mode == GBlendMode::kXor) {
            mode = GBlendMode::kSrcOut;
        }
    }
    
    if (!shader && srcA == 0) {
        if (mode == GBlendMode::kDstOver) {
          return; 
        }
        else if (mode == GBlendMode::kSrcIn) {
            mode = GBlendMode::kSrc;
        }
        else if (mode == GBlendMode::kDstIn) {
            mode = GBlendMode::kSrc;
        }
        else if (mode == GBlendMode::kSrcOut) {
            mode = GBlendMode::kSrc;
        }
        else if (mode == GBlendMode::kDstOut) {
            return;
        }
        else if (mode == GBlendMode::kSrcATop) {
            return;
        }
        else if (mode == GBlendMode::kDstATop) {
            mode = GBlendMode::kSrc;
        }
        else if (mode == GBlendMode::kXor) {
            return;
        }
    }

    GPoint pointsMap[count];
    CTM.mapPoints(pointsMap, points, count);
    vector<Edge> allEdges = createEdges(fDevice, count, pointsMap);

    if (allEdges.size() < 2) {
      return;
    }
      
    sort(allEdges.begin(), allEdges.end(), compare);


    int globalBot = 0;

    for (Edge e : allEdges) {
        if (e.bottom > globalBot) {
            globalBot = e.bottom;
        }
    }

    Edge L = allEdges[0];
    Edge R = allEdges[1];
    
    int cur = 2;

    GPixel *row_addr = nullptr;

    for (int y = max(L.top, 0); y < globalBot; y++) {
        int x1 = GRoundToInt((y +0.5f) *L.m + L.b);
        int x2 = GRoundToInt((y +0.5f) *R.m + R.b);
        int N = x2 - x1;
        if (x1 >= 0 && x1 < fDevice.width() && y >= 0 && y < fDevice.height()) {
          row_addr = fDevice.getAddr(x1, y);
          blitrow(row_addr, N, x1, y, pixel_color, mode, shader);
        } 
        
        if (cur < allEdges.size()) {
            if (!(L.top <= y && y < L.bottom)){
                L = allEdges[cur];
                cur = cur + 1;
            }
        }
        
        if (cur < allEdges.size()) {
            if (!(R.top <= y && y < R.bottom)){
                R = allEdges[cur];
                cur = cur + 1;
            }
        }
    }
    
}


void MyCanvas::drawRect(const GRect& rect, const GPaint& color) {
        GPoint points[4] = {
            { rect.left, rect.bottom },
            { rect.left, rect.top },
            { rect.right, rect.top },
            { rect.right, rect.bottom },
        };
        
        drawConvexPolygon(points, 4, color);
}

std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    return std::unique_ptr<GCanvas>(new MyCanvas(device));
}

std::string GDrawSomething(GCanvas* canvas, GISize dim) {
    // as fancy as you like
    // ...
    // canvas->clear(...);
    /*
    GBitmap bm;
    bm.readFromFile("cr.png");

    GPaint paint(GCreateBitmapShader(bm, GMatrix::Scale(0.5,1)));

    float cx = bm.width();
    float cy = bm.height();
    GPoint pts[] = {
    {0,0}, { cx, 0 }, {cx,cy}, { 0, cy} 
    };
    canvas-> drawConvexPolygon(pts, 4, paint);

    const GPoint pts1[] = {{-30, -30}, {10, 0}, {0, 10}};
    GPathBuilder bu;
    bu.addPolygon(pts1, 3);
    canvas->clear({1,0,1,1});
    canvas->fillRect({1, 1, 8, 8}, {1,1,1,1});
    canvas->fillRect({6, 6, 60, 60}, {1,1,0,1});
    canvas->fillRect({61, 61, 150, 150}, {0,1,1,1});
    */
    GPathBuilder bu;
    bu.reset();
    bu.addCircle({80, 80}, 10, GPathDirection::kCW);
    GRect r = GRect::XYWH(10, 15, 100, 75);
    const GColor colors[] = {
        {1,0,1,1}, {0,1,1,1}, {0,1,1,1},
    };
    auto sh = GCreateLinearGradient({r.left, r.top}, {r.right, r.bottom},
                                    colors, 3);
    GPaint paint(sh);
    canvas->drawRect({-100, -100, 600, 600}, paint); 

    GBitmap bm;
    bm.readFromFile("apps/spock.png");
    float w = bm.width();
    float h = bm.height();

    const GPoint pts[] = {
        { 10, 10 }, { 400, 100 }, { 250, 400 },
    };
    const GPoint tex[] = {
        { 0, 0 }, { w, 0 }, { 0, h },
    };
    const int indices[] = {
        0, 1, 2,
    };
    canvas->drawMesh(pts, nullptr, tex, 1, indices, paint);
    return "Colors";
}

GMatrix::GMatrix() {
  fMat[0] = 1.f;    fMat[2] = 0;    fMat[4] = 0;
  fMat[1] = 0;    fMat[3] = 1.f;    fMat[5] = 0;
}

GMatrix GMatrix::Translate(float tx, float ty) {
  return GMatrix(1, 0, tx,
                 0, 1, ty);
}

GMatrix GMatrix::Scale(float sx, float sy) {
  return GMatrix(sx, 0, 0,
                 0, sy, 0);
}

GMatrix GMatrix::Rotate(float radians) {
  return GMatrix(cos(radians), -sin(radians), 0,
                 sin(radians), cos(radians), 0);
}

GMatrix GMatrix::Concat(const GMatrix& a, const GMatrix& b) {

  float aVal = a.fMat[0] * b.fMat[0] + a.fMat[2] * b.fMat[1];
  float cVal = a.fMat[0] * b.fMat[2] + a.fMat[2] * b.fMat[3];
  float eVal = a.fMat[0] * b.fMat[4] + a.fMat[2] * b.fMat[5] + a.fMat[4];
  float bVal = a.fMat[1] * b.fMat[0] + a.fMat[3] * b.fMat[1];
  float dVal = a.fMat[1] * b.fMat[2] + a.fMat[3] * b.fMat[3];
  float fVal = a.fMat[1] * b.fMat[4] + a.fMat[3] * b.fMat[5] + a.fMat[5];
  return GMatrix(aVal, cVal, eVal, bVal, dVal, fVal);
}


nonstd::optional<GMatrix> GMatrix::invert() const{
  float aVal = fMat[0]; float cVal = fMat[2]; float eVal = fMat[4];
  float bVal = fMat[1]; float dVal = fMat[3]; float fVal = fMat[5];

  float determinant = (aVal * dVal - cVal * bVal);
  if (determinant == 0) {
      return {};
  }

  float invd = 1 / determinant;

  GMatrix inverted = GMatrix(dVal * invd, -1 * cVal * invd, (cVal * fVal - eVal * dVal) * invd,
                             -1 * bVal * invd, aVal * invd, (eVal * bVal - aVal * fVal) * invd);
  return inverted;

}


void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const {
  for (int i=0; i<count; ++i) {
    GPoint g = src[i];

    float x = this->fMat[0] * g.x + this->fMat[2] * g.y + this->fMat[4];
    float y = this->fMat[1] * g.x + this->fMat[3] * g.y + this->fMat[5];

    dst[i] = {x, y};
  }
}


void MyCanvas::save()  {
  stackMatrix.push(CTM);
}

void MyCanvas::restore()  {
  CTM = stackMatrix.top();
  stackMatrix.pop();
}

void MyCanvas::concat(const GMatrix& matrix)  {
  GMatrix creatM = GMatrix::Concat(CTM, matrix);
  CTM = creatM;
}

int numQuad(GPoint pts[3]) {
  GPoint E = (pts[0] - 2 * pts[1] + pts[2]) * 0.25f;
  float mg = sqrt(E.x * E.x + E.y * E.y);

  return (int) ceil(sqrt(mg * 4)); 
}

int numCub(GPoint pts[4]) {
  GPoint E0 = pts[0] - 2 * pts[1] + pts[2];       
  GPoint E1 = pts[1] - 2 * pts[2] + pts[3];         

  float Ex = std::max(abs(E0.x), abs(E1.x)); 
  float Ey = std::max(abs(E0.y), abs(E1.y));
  
  float mg = sqrt(Ex * Ex + Ey * Ey);
  return (int) ceil(sqrt(mg * 3)); 
}

std::vector<Edge> createEdgesPath(const GPath& path, int width, int height) {
  std::vector<Edge> edges;
  GPoint pts[GPath::kMaxNextPoints];
  GPath::Edger edger(path);
  nonstd::optional<GPathVerb> v;

  while ((v = edger.next(pts))) {
      std::vector<Edge> new_edges;
      int seg_num;
      GPoint p0, p1;
      float t;
      GPathVerb g = *v;
      GPoint A;
      GPoint B;
      GPoint C;
      GPoint D;
      switch (g) {
          case GPathVerb::kLine:
              new_edges = clipEdges(height, width, pts[0], pts[1]);

              if (new_edges.size() > 0) {
                edges.insert(edges.end(), new_edges.begin(), new_edges.end());
              }
                              
              break;
          
          case GPathVerb::kCubic:
            seg_num = numCub(pts);
            p0 = pts[0];
            A = pts[0];
            B = pts[1];
            C = pts[2];
            D = pts[3];

            for (int i=1; i < seg_num; i++) {
              t = i * (1.f/seg_num);

              p1 = { getCPt1(A.x, B.x, C.x, D.x, t), getCPt1(A.y, B.y, C.y, D.y, t)};

              new_edges = clipEdges(height, width, p0, p1);

              if (new_edges.size() > 0) {
                  edges.insert(edges.end(), new_edges.begin(), new_edges.end());
              }

              p0 = p1;
            } 

            new_edges = clipEdges(height, width, p0, pts[3]);
            if (new_edges.size() > 0) {
              edges.insert(edges.end(), new_edges.begin(), new_edges.end());
            }

            break;

          case GPathVerb::kQuad:
            seg_num = numQuad(pts);
            p0 = pts[0];

            for (int i=1; i < seg_num; i++) {
              t = i * (1.f/seg_num);
              A = pts[0];
              B = pts[1];
              C = pts[2];

              p1 = { getQPt1(A.x, B.x, C.x, t), getQPt1(A.y, B.y, C.y, t)};

              new_edges = clipEdges(height, width, p0, p1);

              if (new_edges.size() > 0) {
                  edges.insert(edges.end(), new_edges.begin(), new_edges.end());
              }

              p0 = p1;
            } 

            new_edges = clipEdges(height, width, p0, pts[2]);
            if (new_edges.size() > 0) {
              edges.insert(edges.end(), new_edges.begin(), new_edges.end());
            }
            
            break;
          
          default:
           break;

      }
  }
  return edges;
}


bool isValidEdge(Edge e, int y) {
  return e.top <= y && e.bottom > y;
}


void MyCanvas::drawPath(const GPath& path, const GPaint& paint) {
  auto pathMap = path.transform(CTM);

  vector<Edge> allPathEdges = createEdgesPath(*pathMap, fDevice.width(), fDevice.height());


  if (allPathEdges.size() < 2) {
    return;
  }

  sort(allPathEdges.begin(), allPathEdges.end(), compare);

  int L;
  int R;
  int yt = allPathEdges[0].top;
  GBlendMode mode = paint.getBlendMode();
  GPixel pixel_color = ToGPixel(paint.getColor());
  GShader* shader = paint.peekShader();
  GPixel *row_addr = nullptr;

  int globalBot1 = 0;

  for (Edge e : allPathEdges) {
    if (e.bottom > globalBot1) {
      globalBot1 = e.bottom;
    }
  }

  for (int y = yt; y < globalBot1; y++) {
    int w = 0;
    int i = 0;
    while (i<allPathEdges.size() && isValidEdge(allPathEdges[i], y)) {
      if (w == 0) {
        L = GRoundToInt(eval(allPathEdges[i],y));
      }

      assert(allPathEdges[i].wind == 1 || allPathEdges[i].wind == -1);

      w += allPathEdges[i].wind;

      if (w==0) {
        R = GRoundToInt(eval(allPathEdges[i], y));
        int N = R - L;
        if (N > 0) {
        row_addr = fDevice.getAddr(L, y);
        blitrow(row_addr, N, L, y, pixel_color, mode, shader);
        }
      }

      if (!isValidEdge(allPathEdges[i], y+1)) {
        allPathEdges.erase(allPathEdges.begin() + i);
      } else {
        i++;
      }
    }
    
    assert(w == 0);

    while (i < allPathEdges.size() && isValidEdge(allPathEdges[i], y+1)) {
      i++;
    }
    int cy = y + 1;
    std::sort(allPathEdges.begin(), allPathEdges.begin() + i, [cy](Edge a, Edge b) {
      return eval(a, cy) < eval(b, cy); 
    });
  }
}




class BitmapShader : public GShader {
  public:
    BitmapShader(const GBitmap& bitmap, const GMatrix& localMatrix, GTileMode mode)
      : fDevice(bitmap), 
        fLocalMatrix(localMatrix),
        fTileMode(mode) {}
    
    bool isOpaque() override {
      return fDevice.isOpaque();
    }

    bool setContext(const GMatrix& ctm) override {
      if (auto inv = (ctm * fLocalMatrix).invert()) {
        fInverse = *inv;
        return true;
      }
      return false;
    }

    void shadeRow(int x, int y, int count, GPixel row[]) override {
      GPoint p = fInverse * GPoint({ x+0.5f, y+0.5f });

      for (int i = 0; i < count; ++i) {
        int newX;
        int newY; 
        
        switch (fTileMode) {
          case GTileMode::kMirror:
            newX = clamped(mirror(p.x, fDevice.width()), fDevice.width());
            newY = clamped(mirror(p.y, fDevice.height()), fDevice.height());
            break;

          case GTileMode::kRepeat:
            newX = repeat(p.x, fDevice.width());
            newY = repeat(p.y, fDevice.height());
            break;

          case GTileMode::kClamp:
            newX = clamped(p.x, fDevice.width());
            newY = clamped(p.y, fDevice.height());
            break;
        }          
        row[i] = *fDevice.getAddr(newX, newY);

        p.x += fInverse[0];
        p.y += fInverse[1]; 
      }
    }
  private:
    GBitmap fDevice;
    GMatrix fLocalMatrix;
    GMatrix fInverse;
    GTileMode fTileMode;


    int repeat(float point, int dim) {
      float end = point / dim;
      return int((end - GFloorToInt(end)) * dim);
    }

    int mirror(float point, int dim) {
      float end = point / dim;
      int d = GFloorToInt(abs(end));
      if ((d % 2) == 0) {
        return int((end - GFloorToInt(end)) * dim);
      } else {
        return int((GCeilToInt(end) - end) * dim);
      }
    }

    int clamped(float x, int end) {
      return GFloorToInt(x < 0 ? 0: x > end-1 ? end-1: x);
    }

};




std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bitmap, const GMatrix& mat, GTileMode mode) {
  if (!bitmap.pixels()) {
    return nullptr;
  }

  return std::shared_ptr<GShader>(new BitmapShader(bitmap, mat, mode));
}


class LinearGradientShader: public GShader {
public:
	LinearGradientShader(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode) :  fcount(count), fColors(colors, colors+count), fTileMode(mode) {
		for (int i = 0; i < count-1; i++){
			colorsList.push_back(colors[i+1]-colors[i]);
		}

		fCTM = GMatrix();

		float dX = p1.x - p0.x;
		float dY = p1.y - p0.y;

		fLocalMatrix = GMatrix(dX, -dY, p0.x, dY, dX, p0.y);
	}

	bool isOpaque() override {
		for (int i = 0; i < fcount; i++) {
			if (fColors[i].a != 1.0) {
				return false;
			}
		}
		return true;
	}

  bool setContext(const GMatrix& ctm) override {
    if (auto inv = (ctm * fLocalMatrix).invert()) {
      fCTM = *inv;
      return true;
    }
    return false;
  }

	float clamp(GPoint point) {
		float px = point.x;
		if (point.x < 0.0) {
			px = 0;
		}
		if (point.x > 1.0) {
			px = 1;
		}
		return px;
	}

	void shadeRow(int x, int y, int count, GPixel row[]) override {
    GPoint point = { x + 0.5f, y + 0.5f };
    GPoint p = fCTM * point;

    if (fcount == 1) {
      for (int i = 0; i < count; ++i) {
        row[i] = ToGPixel(fColors[0]);
      }
    } else {
      for (int i = 0; i < count; ++i) {
        float x;

        switch (fTileMode) {
          case GTileMode::kMirror:
            x = mirror(p.x);
            break;

          case GTileMode::kRepeat:
            x = repeat(p.x);
            break;

          case GTileMode::kClamp:
            x = clamp(p.x);
            break;
        }
        
        int j = GFloorToInt(x);
        float t = x - j;

        GColor c;
        if (t == 0) {
          c = fColors[j];
        } else {
          c = fColors[j] + (t * colorsList[j]);
        }
        
        row[i] = ToGPixel(c);
        
        p.x += fCTM[0];
      }
    }
    
  }


private:
	int fcount;
	GMatrix fLocalMatrix;
	GMatrix fCTM;
  std::vector<GColor> fColors;
	std::vector<GColor> colorsList;
  GTileMode fTileMode;


  float repeat(float point) {
    return (point - GFloorToInt(point)) * (fcount-1);
  }

  float mirror(float point) {
    int g = GFloorToInt(abs(point));

    if ((g % 2) == 0) {
      return (point - GFloorToInt(point)) * (fcount-1);
    } else {
      return (GCeilToInt(point) - point) * (fcount-1);
    }
  }

  float clamp(float point) {
    return GPinToUnit(point) * (fcount-1);
  }
};


std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode) {
  if (count < 1) {
    return nullptr;
  }

  return std::shared_ptr<GShader>(new LinearGradientShader(p0, p1, colors, count, mode));
}



  class TriangleShader : public GShader {
  public:
    TriangleShader(const GPoint pts[3], const GColor colors[])
      : fColors(colors, colors+3) {

      GPoint U = pts[1] - pts[0];
      GPoint V = pts[2] - pts[0];

      fBoxMatrix =   { U.x, V.x, pts[0].x,
                      U.y, V.y, pts[0].y };

                  
      fColorDiff1 = colors[1] - colors[0];
      fColorDiff2 = colors[2] - colors[0];
    }

    bool isOpaque() override {
      return fColors[0].a == 1.0 && fColors[1].a == 1.0 && fColors[2].a == 1.0;
    }

    bool setContext(const GMatrix& ctm) override {
      if (auto inv = (ctm * fBoxMatrix).invert()) {
        fInverseCTM = *inv;
        return true;
      }
      return false;
    }
    void shadeRow(int x, int y, int count, GPixel row[]) override {
      GPoint point = { x+0.5f, y+0.5f };
      GPoint p = fInverseCTM * point;

      GColor color = p.x * fColorDiff1 + p.y * fColorDiff2 + fColors[0];
      GColor colorChg = fInverseCTM[0] * fColorDiff1 + fInverseCTM[1] * fColorDiff2;
    
      for (int i = 0; i < count; ++i) {
    
        row[i] = clrConversion(color);
        
        p.x += fInverseCTM[0];
        color += colorChg;
      }
    }

    private:

    GPixel clrConversion(GColor color) {
    float a = std::min(std::max(color.a, 0.0f), 1.0f);
    float r = std::min(std::max(color.r * a, 0.0f), 1.0f);
    float g = std::min(std::max(color.g * a, 0.0f), 1.0f);
    float b = std::min(std::max(color.b * a, 0.0f), 1.0f); 

    return GPixel_PackARGB(
        (int)(a * 255.0f + 0.5f),
        (int)(r * 255.0f + 0.5f),
        (int)(g * 255.0f + 0.5f),
        (int)(b * 255.0f + 0.5f)
        );
    }

    std::vector<GColor> fColors;
    GMatrix fInverseCTM;
    GMatrix fBoxMatrix;
    GColor fColorDiff1;
    GColor fColorDiff2;
    int fNumColors;
};

std::shared_ptr<GShader> GCreateTriangleShader(const GPoint pts[3], const GColor colors[]) {
  return std::shared_ptr<GShader>(new TriangleShader(pts, colors));
}


class AdjShader : public GShader {
  public:
      AdjShader(GShader* shader, const GMatrix& bonus)
          : fCurShader(shader), fTransformed(bonus) {}

      bool isOpaque() override { return fCurShader->isOpaque(); }

      bool setContext(const GMatrix& ctm) override {
          return fCurShader->setContext(ctm * fTransformed);
      }
      
      void shadeRow(int x, int y, int count, GPixel row[]) override {
          fCurShader->shadeRow(x, y, count, row);
      }
      
  private:
      GShader* fCurShader;
      GMatrix  fTransformed;
  };

std::shared_ptr<GShader> GCreateAdjShader(GShader* shader, const GMatrix& extraTransform) {
  return std::shared_ptr<GShader>(new AdjShader(shader, extraTransform));
}

class ComboShader : public GShader {
  public:
      ComboShader(GShader* shader0, GShader* shader1)
          : fShader0(shader0), fShader1(shader1) {}

      bool isOpaque() override { return fShader0->isOpaque() && fShader1->isOpaque(); }

      bool setContext(const GMatrix& ctm) override {
          return fShader0->setContext(ctm) && fShader1->setContext(ctm);
      }
      
      void shadeRow(int x, int y, int count, GPixel row[]) override {
          GPixel row0[count];
          GPixel row1[count];

          fShader0->shadeRow(x, y, count, row0);
          fShader1->shadeRow(x, y, count, row1);

          for (int i=0; i<count; i++) {
            row[i] = PixelsMp(row0[i], row1[i]);
          }
      }
      
  private:
      GShader* fShader0;
      GShader* fShader1;

      GPixel PixelsMp(GPixel pixel0, GPixel pixel1) {
        int a = GRoundToInt((GPixel_GetA(pixel0) * GPixel_GetA(pixel1)) / 255.0f);
        int r = GRoundToInt((GPixel_GetR(pixel0) * GPixel_GetR(pixel1)) / 255.0f);
        int g = GRoundToInt((GPixel_GetG(pixel0) * GPixel_GetG(pixel1)) / 255.0f);
        int b = GRoundToInt((GPixel_GetB(pixel0) * GPixel_GetB(pixel1)) / 255.0f);
        
        return GPixel_PackARGB(a, r, g, b);
      }
};

std::shared_ptr<GShader> GCreateComboShader(GShader* shader0, GShader* shader1) {
  return std::shared_ptr<GShader>(new ComboShader(shader0, shader1));
}

GMatrix baseCpt(const GPoint pts[3]) {
  return { pts[1].x - pts[0].x, pts[2].x - pts[0].x, pts[0].x,
            pts[1].y - pts[0].y, pts[2].y - pts[0].y, pts[0].y };
}

void MyCanvas::drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint& paint) {

  int v = 0;
  GShader* shader = paint.peekShader();

  for (int i = 0; i < count; ++i) {
    const GPoint ptsTrig[3] = {verts[indices[v]], verts[indices[v+1]], verts[indices[v+2]]};

    if (texs != nullptr && colors != nullptr) {
      const GColor clrsTrig[3] = { colors[indices[v]], colors[indices[v+1]], colors[indices[v+2]] };
      const GPoint ptsTex[3] = { texs[indices[v]], texs[indices[v+1]], texs[indices[v+2]] };
      GMatrix H;
      GMatrix J;
      GMatrix CT;
      H = baseCpt(ptsTrig);
      J = baseCpt(ptsTex);
      auto W = J.invert();
      CT = *W;

      auto sb = GCreateComboShader(new TriangleShader(ptsTrig, clrsTrig), new AdjShader(shader, H * CT));
      GPaint paint(sb);

      drawConvexPolygon(ptsTrig, 3, paint);


    } else if (colors != nullptr) {
      const GColor clrsTrig[3] = { colors[indices[v]], colors[indices[v+1]], colors[indices[v+2]] };
      auto sf = GCreateTriangleShader(ptsTrig, clrsTrig);
      GPaint paint(sf);

      drawConvexPolygon(ptsTrig, 3, paint);

    } else if (texs != nullptr & shader != nullptr) {
      const GPoint ptsTex[3] = { texs[indices[v]], texs[indices[v+1]], texs[indices[v+2]] };

      GMatrix E;
      GMatrix F;
      GMatrix CG;
      E = baseCpt(ptsTrig);
      F = baseCpt(ptsTex);
      auto B = F.invert();
      CG = *B;
      auto sg = GCreateAdjShader(shader, E * CG);
      GPaint paint(sg);

      drawConvexPolygon(ptsTrig, 3, paint);
    }

  v += 3;
  }

}

GPoint getPtChopped(const GPoint pts[4], float u, float v) {
    return (1 - v) * ((1 - u) * pts[0] + u * pts[1]) +  v * ((1 - u) * pts[3] + u * pts[2]);
}

GColor getClrChopped(const GColor colors[4], float u, float v) {
    float a = (1 - v) * ((1 - u) * colors[0].a + u * colors[1].a) + v * ((1 - u) * colors[3].a + u * colors[2].a);
    float r = (1 - v) * ((1 - u) * colors[0].r + u * colors[1].r) + v * ((1 - u) * colors[3].r + u * colors[2].r);
    float g = (1 - v) * ((1 - u) * colors[0].g + u * colors[1].g) + v * ((1 - u) * colors[3].g + u * colors[2].g);
    float b = (1 - v) * ((1 - u) * colors[0].b + u * colors[1].b) + v * ((1 - u) * colors[3].b + u * colors[2].b);

    return GColor::RGBA(r, g, b, a);
}

void MyCanvas::drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint& paint) {
        GPoint vertsChop[4];
        GColor clrsChop[4];
        GPoint texsChop[4];
    
        for (int u=0; u<=level; u++) {
                float u0 = (float) u / (float) (level + 1.0f);
                float u1 = (float) (u+1) / (float) (level + 1.0f);
            for (int v=0; v<=level; v++) {
                float v0 = (float) v / (float) (level + 1.0f);
                float v1 = (float) (v+1) / (float) (level + 1.0f);
                vertsChop[0] = getPtChopped(verts, u0, v0);
                vertsChop[1] = getPtChopped(verts, u1, v0);
                vertsChop[2] = getPtChopped(verts, u1, v1);
                vertsChop[3] = getPtChopped(verts, u0, v1);

                if (colors != nullptr) {
                    clrsChop[0] = getClrChopped(colors, u0, v0);
                    clrsChop[1] = getClrChopped(colors, u1, v0);
                    clrsChop[2] = getClrChopped(colors, u1, v1);
                    clrsChop[3] = getClrChopped(colors, u0, v1);
                }
                if (texs != nullptr) {
                    texsChop[0] = getPtChopped(texs, u0, v0);
                    texsChop[1] = getPtChopped(texs, u1, v0);
                    texsChop[2] = getPtChopped(texs, u1, v1);
                    texsChop[3] = getPtChopped(texs, u0, v1);
                }
                const int indices[6] = { 0, 1, 3, 1, 2, 3};

                if (colors != nullptr && texs != nullptr) {
                    drawMesh(vertsChop, clrsChop, texsChop, 2, indices, paint);
                } else if (colors != nullptr) {
                    drawMesh(vertsChop, clrsChop, nullptr, 2, indices, paint);
                } else if (texs != nullptr) {
                    drawMesh(vertsChop, nullptr, texsChop, 2, indices, paint);
                } else {
                    drawMesh(vertsChop, nullptr, nullptr, 2, indices, paint);
                }

            }
        }
    }

class LinearPosGradient: public GShader {
public:
	LinearPosGradient(GPoint p0, GPoint p1, const GColor colors[], const float pos[], int count) :  fcount(count), fColors(colors, colors+count), fTileMode(GTileMode::kClamp) {
		for (int i = 0; i < count-1; i++){
			colorsList.push_back(colors[i+1]-colors[i]);
		}

		fCTM = GMatrix();

		float dX = p1.x - p0.x;
		float dY = p1.y - p0.y;

		fLocalMatrix = GMatrix(dX, -dY, p0.x, dY, dX, p0.y);
	}

	bool isOpaque() override {
		for (int i = 0; i < fcount; i++) {
			if (fColors[i].a != 1.0) {
				return false;
			}
		}
		return true;
	}

  bool setContext(const GMatrix& ctm) override {
    if (auto inv = (ctm * fLocalMatrix).invert()) {
      fCTM = *inv;
      return true;
    }
    return false;
  }

	float clamp(GPoint point) {
		float px = point.x;
		if (point.x < 0.0) {
			px = 0;
		}
		if (point.x > 1.0) {
			px = 1;
		}
		return px;
	}

	void shadeRow(int x, int y, int count, GPixel row[]) override {
    GPoint point = { x + 0.5f, y + 0.5f };
    GPoint p = fCTM * point;

    if (fcount == 1) {
      for (int i = 0; i < count; ++i) {
        row[i] = ToGPixel(fColors[0]);
      }
    } else {
      for (int i = 0; i < count; ++i) {
        float x;

        switch (fTileMode) {
          case GTileMode::kMirror:
            x = mirror(p.x);
            break;

          case GTileMode::kRepeat:
            x = repeat(p.x);
            break;

          case GTileMode::kClamp:
            x = clamp(p.x);
            break;
        }
        
        int j = GFloorToInt(x);
        float t = x - j;

        GColor c;
        if (t == 0) {
          c = fColors[j];
        } else {
          c = fColors[j] + (t * colorsList[j]);
        }
        
        row[i] = ToGPixel(c);
        
        p.x += fCTM[0];
      }
    }
    
  }


private:
	int fcount;
	GMatrix fLocalMatrix;
	GMatrix fCTM;
  std::vector<GColor> fColors;
	std::vector<GColor> colorsList;
  GTileMode fTileMode;


  float repeat(float point) {
    return (point - GFloorToInt(point)) * (fcount-1);
  }

  float mirror(float point) {
    int g = GFloorToInt(abs(point));

    if ((g % 2) == 0) {
      return (point - GFloorToInt(point)) * (fcount-1);
    } else {
      return (GCeilToInt(point) - point) * (fcount-1);
    }
  }

  float clamp(float point) {
    return GPinToUnit(point) * (fcount-1);
  }
};


class TheFinal : public GFinal {
public:
    TheFinal() {}

    std::shared_ptr<GShader> createLinearPosGradient(GPoint p0, GPoint p1,
                                                             const GColor colors[],
                                                             const float pos[],
                                                             int count) {
        return std::shared_ptr<GShader>(new LinearPosGradient(p0, p1, colors, pos, count));
    }

    std::shared_ptr<GPath> strokePolygon(const GPoint points[], int count, float width, bool isClosed) {
 
      GPathBuilder strokePath;

      for (int i=0; i<count-1; i++) {
        GPoint p0 = points[i];
        GPoint p1 = points[i+1];
        if (p0.x > p1.x) { std::swap(p0, p1); }

        float delta_x = p1.x - p0.x;
        float delta_y = p1.y - p0.y;

        float dist = width/2.0;

        GPoint check = {-delta_y,delta_x};
        check = check * (1.0/sqrt(check.x * check.x + check.y * check.y));

        strokePath.moveTo(p0 + check * dist);

        strokePath.lineTo(p0 - check * dist);
        strokePath.lineTo(p1 - check * dist);
        strokePath.lineTo(p1 + check * dist);

        strokePath.addCircle(p0, dist);
        strokePath.addCircle(p1, dist);
        }

      if (isClosed) {
        GPoint p0 = points[0];
        GPoint p1 = points[count-1];
        if (p0.x > p1.x) { std::swap(p0, p1); }

        float delta_x2 = p1.x - p0.x;
        float delta_y2 = p1.y - p0.y;


        float dist2 = width/2.0;

        GPoint check2 = {-delta_y2,delta_x2};
        check2 = check2 * (1.0/sqrt(check2.x * check2.x + check2.y * check2.y));

        strokePath.moveTo(p0 + check2 * dist2);

        strokePath.lineTo(p0 - check2 * dist2);
        strokePath.lineTo(p1 - check2 * dist2);
        strokePath.lineTo(p1 + check2 * dist2);

        strokePath.addCircle(p0, dist2);
        strokePath.addCircle(p1, dist2);
      }

      return strokePath.detach();
    }
};


std::unique_ptr<GFinal> GCreateFinal() {
  return std::make_unique<TheFinal>();
}
