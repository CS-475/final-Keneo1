#include "include/GMatrix.h"
#include "include/GPath.h"
#include "include/GPathBuilder.h"
#include <iostream>
using namespace std;


float getQPt(float A, float B, float C, float t) {
  float qpt = A + t * (B - A) + t * ((B + t * (C - B)) - (A + t * (B - A)));

  return qpt;
}

float getCPt(float A, float B, float C, float D, float t) {
  float quad1 = getQPt(A, B, C, t);
  float quad2 = getQPt(B, C, D, t);
  float cubic1 = quad1 + t * (quad2 - quad1);

  return cubic1;
}


void GPathBuilder::addPolygon(const GPoint pts[], int count) {
    assert(count >= 2);
    this->moveTo(pts[0]);

    for (int i = 1; i < count; i++) {
        this->lineTo(pts[i]);
    }

}

void GPathBuilder::addRect(const GRect& rect, GPathDirection dir) {

    this->moveTo({rect.left, rect.top});

    if (dir == GPathDirection::kCW) {
        this->lineTo({rect.right, rect.top});
        this->lineTo({rect.right, rect.bottom});
        this->lineTo({rect.left, rect.bottom});
  } else {
        this->lineTo({rect.left, rect.bottom});
        this->lineTo({rect.right, rect.bottom});
        this->lineTo({rect.right, rect.top});
  }
}


GRect GPath::bounds() const {  
  int count = this->fPts.size();

  if (count == 0) {
    return GRect::LTRB(0,0,0,0);
  }  

  if (count == 1) {
    return GRect::XYWH(fPts[0].x, fPts[0].y, 0, 0);
  }

  float minX = fPts[0].x;
  float maxX = fPts[0].x;
  float minY = fPts[0].y;
  float maxY = fPts[0].y;

  GPoint pts[GPath::kMaxNextPoints];
  GPath::Edger iter(*this);
  nonstd::optional<GPathVerb> v;
  while ((v = iter.next(pts)) ) {
      GPoint A, B, C, D;
      float t0;
      float t1;
      float ax;
      float bx;
      float cx;
      float ay;
      float by;
      float cy;

      float px0 = minX;
      float py0 = minY;
      float px1 = minX; 
      float py1 = minY;

      GPathVerb g = *v;
      switch (g) {
        case GPathVerb::kLine:
          minX = min({ minX, pts[0].x, pts[1].x });
          maxX = max({ maxX, pts[0].x, pts[1].x } );
          minY = min({ minY, pts[0].y, pts[1].y });
          maxY = max({ maxY, pts[0].y, pts[1].y });

          break;

        case GPathVerb::kCubic:
          A = pts[0];
          B = pts[1];
          C = pts[2];
          D = pts[3];

          ax = -A.x + 3 * B.x - 3 * C.x + D.x;
          bx = 2 * A.x - 4 * B.x + 2 * C.x;
          cx = -A.x + B.x;
          ay = -A.y + 3 * B.y - 3 * C.y + D.y;
          by = 2 * A.y - 4 * B.y + 2 * C.y;
          cy = -A.y + B.y;

          if (ax == 0) {
            t0 = -cx / bx;
            t1 = -1;
          } else {
            t0 = (-bx + sqrt(bx*bx - 4 * ax * cx)) / (2 * ax);
            t1 = (-bx - sqrt(bx*bx - 4 * ax * cx)) / (2 * ax);
          }
          
          if(t1 >= 0 && t1 <= 1) {
            px1 = getCPt(A.x, B.x, C.x, D.x, t1);
            py1 = getCPt(A.y, B.y, C.y, D.y, t1);
          }

          if(t0 >= 0 && t0 <= 1) {
            px0 = getCPt(A.x, B.x, C.x, D.x, t0);
            py0 = getCPt(A.y, B.y, C.y, D.y, t0);
          }

          minX = min({ minX, px0, px1, A.x, D.x });
          maxX = max({ maxX, px0, px1, A.x, D.x });
          minY = min({ minY, py0, py1, A.y, D.y });
          maxY = max({ maxY, py0, py1, A.y, D.y });

          if (ay == 0) {
            t0 = -cy / by;
            t1 = -1;
          } else {
            t0 = (-by + sqrt(by*by - 4 * ay * cy)) / (2 * ay);
            t1 = (-by - sqrt(by*by - 4 * ay * cy)) / (2 * ay);
          }

          if(t1 >= 0 && t1 <= 1) {
            px1 = getCPt(A.x, B.x, C.x, D.x, t1);
            py1 = getCPt(A.y, B.y, C.y, D.y, t1);
          }

          if(t0 >= 0 && t0 <= 1) {
            px0 = getCPt(A.x, B.x, C.x, D.x, t0);
            py0 = getCPt(A.y, B.y, C.y, D.y, t0);
          }

          minX = min({ minX, px0, px1, A.x, D.x });
          maxX = max({ maxX, px0, px1, A.x, D.x });
          minY = min({ minY, py0, py1, A.y, D.y });
          maxY = max({ maxY, py0, py1, A.y, D.y });

          break;

        case GPathVerb::kQuad:
          A = pts[0];
          B = pts[1];
          C = pts[2];
          t0 = (A.x - B.x) / (A.x - 2 * B.x + C.x);
          t1 = (A.y - B.y) / (A.y - 2 * B.y + C.y);

          if(t0 >= 0 && t0 <= 1) {
            px0 = getQPt(A.x, B.x, C.x, t0);
            py0 = getQPt(A.y, B.y, C.y, t0);
          }

          if(t1 >= 0 && t1 <= 1) {
            px1 = getQPt(A.x, B.x, C.x, t1);
            py1 = getQPt(A.y, B.y, C.y, t1);
          }

          minX = min({ minX, px0, px1, A.x, C.x });
          maxX = max({ maxX, px0, px1, A.x, C.x });
          minY = min({ minY, py0, py1, A.y, C.y });
          maxY = max({ maxY, py0, py1, A.y, C.y });

          break;

        default:
          break;
      }
  }
  return GRect{minX,minY,maxX,maxY};
}



void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection dir) {
    float x = center.x;
    float y = center.y;

    float c = 0.551915f;
    float distance = c * radius;
    this->moveTo(x + radius, y);

  if (dir == GPathDirection::kCW) {
    this->cubicTo({ x + radius, y + distance }, { x + distance, y + radius }, { x, y + radius }); 
    this->cubicTo({ x - distance, y + radius }, { x - radius, y + distance }, { x - radius, y }); 
    this->cubicTo({ x - radius, y - distance }, { x - distance, y - radius }, { x, y - radius});  
    this->cubicTo({ x + distance, y - radius }, { x + radius, y - distance }, { x + radius, y });
  } else {
    this->cubicTo({ x + radius, y - distance }, { x + distance, y - radius }, { x, y - radius});
    this->cubicTo({ x - distance, y - radius }, { x - radius, y - distance }, { x - radius, y });
    this->cubicTo({ x - radius, y + distance }, { x - distance, y + radius }, { x, y + radius });
    this->cubicTo({ x + distance, y + radius }, { x + radius, y + distance }, { x + radius, y });
  }
}


void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t) {
    dst[0] = src[0];
    dst[1] = (1 - t) * src[0] + t * src[1];
    dst[2] = (1 - t) * (1 - t) * src[0] + 2 * t * (1 - t) * src[1] + t * t * src[2];
    dst[3] = (1 - t) * src[1] + t * src[2];
    dst[4] = src[2];
}

void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t) {
    dst[0] = src[0];
    dst[1] = (1 - t) * src[0] + t * src[1];
    dst[2] = (1 - t) * ((1 - t) * src[0] + t * src[1]) + t * ((1 - t) * src[1] + t * src[2]);
    dst[3] = (1 - t) * (1 - t) * (1 - t) * src[0] + 3 * t * (1 - t) * (1 - t) * src[1] + 3 * t * t * (1 - t) * src[2] + t * t * t * src[3];
    dst[4] = (1 - t) * ((1 - t) * src[1] + t * src[2]) + t * ((1 - t) * src[2] + t * src[3]);
    dst[5] = (1 - t) * src[2] + t * src[3];
    dst[6] = src[3];
}