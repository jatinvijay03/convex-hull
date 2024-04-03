import 'dart:math';
import 'Custom_Point.dart';

animatedPointsjm jmfin = animatedPointsjm([], []); /// Initialize animation data

/// Function to determine the orientation of three points (p, q, r). It uses the slope format and some mathematical simplification in order to generate an equation to find the direction of another point (in this case r) with respect to 2 given points (in this case p and q) 
/// Returns:
///   0 if points are collinear
///   1 if points are in a clockwise orientation
///   2 if points are in a counterclockwise orientation
int orientation(Custom_Point p, Custom_Point q, Custom_Point r) {
  double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

  if (val == 0) return 0; /// Collinear
  return (val > 0) ? 1 : 2; /// Clockwise or counterclockwise
}

/// Function to calculate the squared distance between two points
double distanceSquared(Custom_Point p, Custom_Point q) {
  return (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y);
}

/// Function to check if a point is inside a polygon. It uses a ray tracing method, where it draws an infinite ray toward the right from the point. The number of times the ray intersects with the polygon's sides is counted. The point is inside the polygon if the number of intersections is odd. If the number of intersections is even, then the point is outside the polygon.
bool pointInPolygon(Custom_Point point, List<Custom_Point> polygon) {
  int numVertices = polygon.length;
  double x = point.x, y = point.y;
  bool inside = false;
  Custom_Point p1 = polygon[0], p2;

  for (int i = 1; i <= numVertices; i++) {
    p2 = polygon[i % numVertices];
    if (y > min(p1.y, p2.y)) {
      if (y <= max(p1.y, p2.y)) {
        if (x <= max(p1.x, p2.x)) {
          double xIntersection =
              (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
          if (p1.x == p2.x || x <= xIntersection) {
            inside = !inside;
          }
        }
      }
    }
    p1 = p2;
  }
  return inside;
}

/// Function to compute the convex hull of a set of points using Jarvis March algorithm
List<Custom_Point> convexHull(List<Custom_Point> points) {
  int n = points.length;
  if (n < 3) return []; /// Convex hull not possible with less than 3 points

  List<Custom_Point> hull = []; /// Initialize hull

  /// Find the leftmost point
  int l = 0;
  for (int i = 1; i < n; i++) {
    if (points[i].x < points[l].x) {
      l = i;
    } else if (points[i].x == points[l].x && points[i].y < points[l].y) {
      l = i;
    }
  }

  int p = l, q;
  /// Start Jarvis March algorithm
  do {
    hull.add(points[p]); /// Add current point to the hull
    q = (p + 1) % points.length; /// Set next point to the next point in array

    /// Find the next point among the remaining points which have a possibility of being on a more counterclosckwise / clockwise direction, depending upon the need (in our case it is counterclockwise)
    for (int i = 0; i < points.length; i++) {
      if (orientation(points[p], points[i], points[q]) == 2) {
        q = i;
      }
    }

    List<Custom_Point> pointsToRemove = [];
    List<Custom_Point> pointsCopy = List.from(points);

    /// Check if there are enough points to form a convex hull
    if (hull.length >= 3) {
      for (int i = 0; i < pointsCopy.length; i++) {
        if (pointInPolygon(pointsCopy[i], hull)) {
          /// Remove points inside the current convex hull
          pointsToRemove.add(pointsCopy[i]);
        }
      }
      for (int i = 0; i < hull.length; i++) {
        pointsToRemove.add(hull[i]);
      }

      /// Store removed points for animation
      jmfin.removedP.add(pointsToRemove);

      /// Remove points from the original list
      for (int i = 0; i < pointsToRemove.length; i++) {
        pointsCopy.remove(pointsToRemove[i]);
      }
    }
    p = q; /// Set next point as current point
  } while (p != l); /// Repeat until the first point is reached again

  hull.add(hull[0]); /// Close the hull
  return hull;
}

/// Function to compute the convex hull of a set of points using Jarvis March algorithm
animatedPointsjm jarvisMarch(List<Custom_Point> points) {
  convexHull(points); /// Compute convex hull
  return jmfin; /// Return animation data
}