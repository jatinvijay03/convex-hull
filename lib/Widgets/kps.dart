import 'dart:math';
import 'Custom_Point.dart';

/// The Line class represents a line segment defined by two points (p1, p2).
/// It includes a method to calculate the slope of the line, which is essential for comparing lines.
class Line {
  Custom_Point p1;
  Custom_Point p2;

  Line(this.p1, this.p2);

  double slope() {
    if (p1.x == p2.x) {
      return double.infinity;
    }
    return (p2.y - p1.y) / (p2.x - p1.x);
  }

  @override
  String toString() => 'Line from $p1 to $p2';
}

List <List<Custom_Point>> upperBridgePoints = [];
List <List<Custom_Point>> lowerBridgePoints = [];
List <List<Custom_Point>> quadrilateral = [];
List <List<Custom_Point>> removedPoints = [];
List <Custom_Point> cvhull = [];
List <List<List<Custom_Point>>> Ordered = [];


findkps(points){
  upperBridgePoints = [];
  lowerBridgePoints = [];
  quadrilateral = [];
  removedPoints = [];
  cvhull = [];
  Ordered = [];
  kirkPatrick(points);
  animatedPoints fin = animatedPoints(upperBridgePoints,lowerBridgePoints,quadrilateral,removedPoints,cvhull,Ordered);
  return fin;
}

/// Determined if a point is inside a set of points(polygon here)
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

/// Function to remove consecutve duplicate entries from a list
List<Custom_Point> removeConsecutiveDuplicates(List<Custom_Point> points) {
  return points
      .asMap()
      .entries
      .where((entry) => 
        entry.key == 0 || points[entry.key] != points[entry.key - 1])
      .map((entry) => entry.value)
      .toList();
}

/// O(n) function to help determine the median of a list. Helps reduce the overall time complexity of the algorithm as calculation of median is requied in many places.
double medianOfMedians(List<double> list, int i) {
  if (list.length <= 5) {
    list.sort();
    return list[i % list.length]; 
  }
  /// Divide the list into groups of 5, sort each group, and find their medians.
  List<List<double>> groups = [];
  for (int j = 0; j < list.length; j += 5) {
    List<double> group = list.sublist(j, min(j + 5, list.length));
    group.sort();
    groups.add(group);
  }
  /// Recursively find the median of the medians to use as a pivot for partitioning.
  List<double> medians = groups.map((group) => group[group.length ~/ 2]).toList();

  double pivot = medianOfMedians(medians, medians.length ~/ 2);
  /// Partition the original list around the pivot and recursively find the kth smallest element if necessary.
  List<double> less = [];
  List<double> equal = [];
  List<double> greater = [];
  for (double num in list) {
    if (num < pivot) {
      less.add(num);
    } else if (num > pivot) {
      greater.add(num);
    } else {
      equal.add(num);
    }
  }


  if (i < less.length) {
    return medianOfMedians(less, i);  
  } else if (i < less.length + equal.length) {
    return pivot;  
  } else {
    return medianOfMedians(greater, i - less.length - equal.length);  
  }
}

/// Finds the point with the minimum x-coordinate (and maximum y if there's a tie).
/// This point is one of the starting points for constructing the convex hull.
Custom_Point findPumin(List<Custom_Point> points) {
  Custom_Point pumin = points[0];
  for (Custom_Point p in points) {
    if (p.x < pumin.x || (p.x == pumin.x && p.y > pumin.y)) {
      pumin = p;
    }
  }
  return pumin;
}

/// Similar to findPumin but for the maximum x-coordinate.
Custom_Point findPumax(List<Custom_Point> points) {
  Custom_Point pumax = points[0];
  for (Custom_Point p in points) {
    if (p.x > pumax.x || (p.x == pumax.x && p.y > pumax.y)) {
      pumax = p;
    }
  }
  return pumax;
}
/// Finds the point with the minimum x-coordinate but the minimum y-coordinate for ties.
/// This is used for the lower hull computation.
Custom_Point findPlmin(List<Custom_Point> points) {
  Custom_Point plmin = points[0];
  for (Custom_Point p in points) {
    if (p.x < plmin.x || (p.x == plmin.x && p.y < plmin.y)) {
      plmin = p;
    }
  }
  return plmin;
}
/// Similar to findPlmin but for the maximum x-coordinate and minimum y in ties, for the lower hull.
Custom_Point findPlmax(List<Custom_Point> points) {
  Custom_Point plmax = points[0];
  for (Custom_Point p in points) {
    if (p.x > plmax.x || (p.x == plmax.x && p.y < plmax.y)) {
      plmax = p;
    }
  }
  return plmax;
}

/// Computes a line (medianLine) that represents the median x-coordinate of the points.
/// This line is used to divide the point set into two subsets for recursive computation.
Line findMedianLine(List<Custom_Point> points) {
  List<double> xCoordinates = points.map((point) => point.x).toList();
  double medianX = medianOfMedians(xCoordinates, xCoordinates.length ~/ 2);
  Custom_Point medianPoint1 = Custom_Point(medianX, double.negativeInfinity);
  Custom_Point medianPoint2 = Custom_Point(medianX, double.infinity);
  return Line(medianPoint1, medianPoint2);
}

/// Divides points into two lists based on their position relative to the median line.
void splitPoints(List<Custom_Point> points, Line medianLine, List<Custom_Point> tLeft, List<Custom_Point> tRight) {
  double medianX = medianLine.p1.x; 
  for (Custom_Point point in points) {
    if (point.x < medianX) {
      tLeft.add(point);
    } else if (point.x > medianX) {
      tRight.add(point);
    } else {
      tLeft.add(point);
    }
  }
}
/// Finds the upper bridge of the point set, which is a line between two points that divides
/// the set such that all other points lie below it.
/// This function is critical for constructing the upper hull of the convex hull.
Line findUpperBridge(List<Custom_Point> points, Line medianLine) {
   /// Base case: directly return the line between two points.
  if (points.length == 2) {
    return Line(points[0], points[1]);
  }

  List<Line> pairs = [];
  List<Custom_Point> candidates = [];
  for (int i = 0; i < points.length - 1; i += 2) {
    pairs.add(Line(points[i], points[i + 1]));
  }
  if (points.length % 2 != 0) { /// Handle an odd number of points.
    candidates.add(points.last);
  }
   /// Calculate slopes of each pair to find the median slope, which guides the selection of the upper bridge.
  List<double> slopes = [];
  List<Line> noninfslopes = [];
  
  List<Line> smallSlopePairs = [];
  List<Line> equalSlopePairs = [];
  List<Line> largeSlopePairs = [];
  
  for(Line pair in pairs) {
    if(pair.p1.x == pair.p2.x) {
      if(pair.p1.y > pair.p2.y) { /// For vertical lines, the higher point is a candidate.
        candidates.add(pair.p1);
      }
      else {
        candidates.add(pair.p2);
      }
    }
    else {
      double pairSlope = pair.slope();
      slopes.add(pairSlope);
      noninfslopes.add(pair);
    }
  }
  double medianSlope = medianOfMedians(slopes, slopes.length ~/ 2);
  /// Categorize non-vertical pairs into three lists based on their slope relative to the median slope.
  for(Line pair in noninfslopes) {
      double pairSlope = pair.slope();
      if (pairSlope < medianSlope) {
        smallSlopePairs.add(pair);
      } 
      else if (pairSlope == medianSlope) {
        equalSlopePairs.add(pair);
      } 
      else { 
        largeSlopePairs.add(pair);
      }
    }
  
  /// Find the highest point on the median slope line, which will be part of the upper bridge.
  double maxYIntercept = -double.infinity;
  Custom_Point? temp;
  Custom_Point pk = Custom_Point(double.infinity, 0);
  Custom_Point pm = Custom_Point(-double.infinity, 0);

  for(Custom_Point point in points)
  {
    double yIntercept = point.y - medianSlope * point.x;
    if (yIntercept > maxYIntercept) {
      maxYIntercept=yIntercept;
    }
  }

  for(Custom_Point point in points)
  {
    double yIntercept = point.y - medianSlope * point.x;
    if (yIntercept == maxYIntercept) {
      if(point.x<pk.x)
      {
        pk=point;
      }
      if(point.x>pm.x)
      {
        pm=point;
      }
    }
  }

  if (pk == null || pm == null) {
    throw Exception('No supporting line found.');
  }
  /// Recursively find the upper bridge if it's not directly determined by the current candidates.
  double a = medianLine.p1.x;
  if (pk.x <= a && pm.x > a) {
    return Line(pk, pm);
  } else {
    if (pm.x <= a) { 
      /// If all points of the upper bridge are to the left of the median, include points from LARGE and EQUAL slopes, and both points from SMALL slope pairs.
      candidates.addAll(largeSlopePairs.map((pair) => pair.p2));
      candidates.addAll(equalSlopePairs.map((pair) => pair.p2));
      candidates.addAll(smallSlopePairs.expand((pair) => [pair.p1, pair.p2]));
    } 
    else if(pk.x > a) {
      /// If all points of the upper bridge are to the right of the median, include points from SMALL and EQUAL slopes, and both points from LARGE slope pairs.
      candidates.addAll(smallSlopePairs.map((pair) => pair.p1));
      candidates.addAll(equalSlopePairs.map((pair) => pair.p1));
      candidates.addAll(largeSlopePairs.expand((pair) => [pair.p1, pair.p2]));
    }
  }
  candidates.sort((a, b) => a.x.compareTo(b.x));
  return findUpperBridge(candidates, medianLine);  
}
/// Similar to findUpperBridge, but for finding the lower bridge.
/// The logic is mirrored to ensure all points lie above the bridge, crucial for the lower hull construction.
Line findLowerBridge(List <Custom_Point> points, Line medianLine) {
   /// Base case: directly return the line between two points.
  if(points.length == 2){
    return Line(points[0], points[1]);
  }
    
  List<Line> pairs = [];
  List<Custom_Point> candidates = [];
  for (int i = 0; i < points.length - 1; i += 2) {
    pairs.add(Line(points[i], points[i + 1]));
  }
  if (points.length % 2 != 0) { /// Handle an odd number of points.
    candidates.add(points.last);
  }
  /// Calculate slopes for the median slope selection, excluding vertical lines.
  List<double> slopes = [];
  List<Line> noninfslopes = [];
  
  List<Line> smallSlopePairs = [];
  List<Line> equalSlopePairs = [];
  List<Line> largeSlopePairs = [];

  for(Line pair in pairs) {
    if(pair.p1.x == pair.p2.x) {
      /// For vertical lines, the lower point is a candidate for the lower bridge.
      if(pair.p1.y < pair.p2.y) {
        candidates.add(pair.p1);
      }
      else {
        candidates.add(pair.p2);
      }
    }
    else {
      double pairSlope = pair.slope();
      slopes.add(pairSlope);
      noninfslopes.add(pair);
    }
  }
  double medianSlope = medianOfMedians(slopes, slopes.length ~/ 2);
  /// Categorize lines based on their slope relative to the median slope.
  for(Line pair in noninfslopes) {
      double pairSlope = pair.slope();
      if (pairSlope < medianSlope) {
        smallSlopePairs.add(pair);
      } 
      else if (pairSlope == medianSlope) {
        equalSlopePairs.add(pair);
      } 
      else { 
        largeSlopePairs.add(pair);
      }
    }
  
  /// Find the lowest point on the median slope line, which will be part of the lower bridge.
  double minYIntercept = double.infinity;
  Custom_Point pk = Custom_Point(double.infinity, 0);
  Custom_Point pm = Custom_Point(-double.infinity, 0);

  for (Custom_Point point in points) {
    double yIntercept = point.y - medianSlope * point.x;
    if (yIntercept < minYIntercept) {
      minYIntercept = yIntercept;
    }
  }

  for (Custom_Point point in points) {
    double yIntercept = point.y - medianSlope * point.x;
    if (yIntercept == minYIntercept) {
      if(point.x<pk.x)
      {
        pk=point;
      }
      if(point.x>pm.x)
      {
        pm=point;
      }
    }
  }

  if (pk == null || pm == null) {
    throw Exception('No supporting line found.');
  }
  /// Recursively find the lower bridge if not directly determined.
  double a = medianLine.p1.x;
  if (pk.x <= a && pm.x > a) {
    return Line(pk, pm);
  }
  else {
    if (pm.x <=a) {  
      /// If all points of the lower bridge are to the left of the median, include points from SMALL slopes and both points from LARGE slope pairs.
      candidates.addAll(smallSlopePairs.map((pair) => pair.p2));
      candidates.addAll(equalSlopePairs.map((pair) => pair.p2));
      candidates.addAll(largeSlopePairs.expand((pair) => [pair.p1, pair.p2]));
    } 
    else if(pm.x > a) {
      /// If all points of the lower bridge are to the right of the median, include points from LARGE slopes and both points from SMALL slope pairs.
      candidates.addAll(largeSlopePairs.map((pair) => pair.p1));
      candidates.addAll(equalSlopePairs.map((pair) => pair.p1));
      candidates.addAll(smallSlopePairs.expand((pair) => [pair.p1, pair.p2]));
    }
  }
  candidates.sort((a, b) => a.x.compareTo(b.x));
  return findLowerBridge(candidates, medianLine);
}

/// Constructs the upper hull of the convex hull by recursively finding upper bridges and removing points within the hull.
/// This is a divide-and-conquer approach that builds the upper part of the convex hull.
List<Custom_Point> upperHull(Custom_Point pUmin, Custom_Point pUmax, List <Custom_Point> points) {
  /// Base case: the points themselves form part of the upper hull.
  if(points.length <= 2){
    upperBridgePoints.add(points);
    Ordered.add([points,[],[]]);
    return points;
  }
  List<Custom_Point> pointsCopy = List.from(points);
  Line medianLine = findMedianLine(points);
  Line ubridge = findUpperBridge(pointsCopy,medianLine);
  List<Custom_Point> temp = [];
  temp.add(pUmin);

  temp.add(ubridge.p1);
  temp.add(ubridge.p2);

  temp.add(pUmax);
  upperBridgePoints.add(temp);
  List<Custom_Point> pointsToRemove = [];
  List<Custom_Point> uBridgePoints = [];
  uBridgePoints.add(pUmin);
  uBridgePoints.add(ubridge.p1);
  uBridgePoints.add(ubridge.p2);
  uBridgePoints.add(pUmax);
  uBridgePoints.add(pUmin);

  quadrilateral.add(uBridgePoints);
  for (int i = 0; i < pointsCopy.length; i++) {
        if (pointInPolygon(pointsCopy[i], uBridgePoints) && !uBridgePoints.contains(pointsCopy[i])) {
          pointsToRemove.add(pointsCopy[i]);
        }
  }
  removedPoints.add(pointsToRemove);
  for (int i = 0; i < pointsToRemove.length; i++) {
        pointsCopy.remove(pointsToRemove[i]);
  }
  /// Divide the remaining points based on their position relative to the upper bridge for recursive computation.
  List<Custom_Point> tLeft = [];
  List<Custom_Point> tRight = [];
  double slopeULeft = (uBridgePoints[1].y - pUmin.y)/(uBridgePoints[1].x - pUmin.x); /// Slope of the line from pUmin to the left point of the upper bridge.
  double yInterceptULeft = uBridgePoints[1].y - (slopeULeft * uBridgePoints[1].x); /// Y-intercept of the same line.
  for(int i = 0; i < pointsCopy.length; i++) {
    double yIntercept = pointsCopy[i].y - (pointsCopy[i].x * slopeULeft);
    if(yIntercept > yInterceptULeft) {
      tLeft.add(pointsCopy[i]);
    }
  }
  tLeft.add(uBridgePoints[1]);
  double slopeURight = (pUmax.y - uBridgePoints[2].y)/(pUmax.x - uBridgePoints[2].x); /// Slope of the line from pUmax to the right point of the upper bridge.
  double yInterceptURight = uBridgePoints[2].y - (slopeURight * uBridgePoints[2].x); /// Y-intercept of the same line.
  for(int i = 0; i < pointsCopy.length; i++) {
    double yIntercept = pointsCopy[i].y - (pointsCopy[i].x * slopeURight);
    if(yIntercept > yInterceptURight) {
      tRight.add(pointsCopy[i]); /// Points above the line are part of the right set.
    }
  }
  tRight.add(uBridgePoints[2]); /// Include the bridge point in the set for recursion.
  /// Recursively construct the upper hull for the left and right subsets.
  List<Custom_Point> upperHullLeft = upperHull(pUmin,uBridgePoints[1],tLeft);
  List<Custom_Point> upperHullRight = upperHull(uBridgePoints[2],pUmax,tRight);
  List<Custom_Point> upperHullFinal = [];
  upperHullFinal.addAll(uBridgePoints);
  upperHullFinal.addAll(upperHullLeft);
  upperHullFinal.addAll(upperHullRight);
  Ordered.add([temp,uBridgePoints,pointsToRemove]);
  return upperHullFinal; /// Return the complete list of points forming the upper hull.
}

/// Constructs the lower hull of the convex hull by a process similar to upperHull.
/// This function handles the bottom part of the convex hull, ensuring all points lie above the constructed hull.
List<Custom_Point> lowerHull(Custom_Point pLmin, Custom_Point pLmax, List <Custom_Point> points) {
  /// Base case: the points themselves form part of the lower hull.
  if(points.length <= 2){
    lowerBridgePoints.add(points);
    Ordered.add([points,[],[]]);
    return points;
  }
  List<Custom_Point> pointsCopy = List.from(points);
  Line medianLine = findMedianLine(points);
  Line lbridge = findLowerBridge(pointsCopy,medianLine);
  List<Custom_Point> temp = [];
  temp.add(pLmin);
  temp.add(lbridge.p1);
  temp.add(lbridge.p2);
  temp.add(pLmax);
  lowerBridgePoints.add(temp);
  List<Custom_Point> pointsToRemove = [];
  List<Custom_Point> lBridgePoints = [];
  lBridgePoints.add(pLmin);
  lBridgePoints.add(lbridge.p1);
  lBridgePoints.add(lbridge.p2);
  lBridgePoints.add(pLmax);
  lBridgePoints.add(pLmin);
  quadrilateral.add(lBridgePoints);
  for (int i = 0; i < pointsCopy.length; i++) {
        if (pointInPolygon(pointsCopy[i], lBridgePoints) && !lBridgePoints.contains(pointsCopy[i])) {
          pointsToRemove.add(pointsCopy[i]);
        }
  }
  removedPoints.add(pointsToRemove);
  Ordered.add([temp,lBridgePoints,pointsToRemove]);
  for (int i = 0; i < pointsToRemove.length; i++) {
        pointsCopy.remove(pointsToRemove[i]);
  }
  /// Divide the remaining points based on their position relative to the lower bridge for recursive computation.
  List<Custom_Point> tLeft = [];
  List<Custom_Point> tRight = [];
  double slopeLLeft = (lBridgePoints[1].y - pLmin.y)/(lBridgePoints[1].x - pLmin.x); /// Slope of the line from pLmin to the left point of the lower bridge.
  double yInterceptLLeft = lBridgePoints[1].y - (slopeLLeft * lBridgePoints[1].x); /// Y-intercept of the same line.
  for(int i = 0; i < pointsCopy.length; i++) {
    double yIntercept = pointsCopy[i].y - (pointsCopy[i].x * slopeLLeft);
    if(yIntercept < yInterceptLLeft) {
      tLeft.add(pointsCopy[i]);
    }
  }
  tLeft.add(lBridgePoints[1]);
  double slopeLRight = (pLmax.y - lBridgePoints[2].y)/(pLmax.x - lBridgePoints[2].x); /// Slope of the line from pLmax to the right point of the lower bridge.
  double yInterceptLRight = lBridgePoints[2].y - (slopeLRight * lBridgePoints[2].x); /// Y-intercept of the same line.
  for(int i = 0; i < pointsCopy.length; i++) {
    double yIntercept = pointsCopy[i].y - (pointsCopy[i].x * slopeLRight);
    if(yIntercept < yInterceptLRight) {
      tRight.add(pointsCopy[i]);
    }
  }
  tRight.add(lBridgePoints[2]);
  List<Custom_Point> lowerHullLeft = lowerHull(pLmax,lBridgePoints[2],tRight);
  List<Custom_Point> lowerHullRight = lowerHull(lBridgePoints[1],pLmin,tLeft);
  List<Custom_Point> lowerHullFinal = [];
  lowerHullFinal.addAll(lBridgePoints);
  lowerHullFinal.addAll(lowerHullLeft);
  lowerHullFinal.addAll(lowerHullRight);
  
  return lowerHullFinal;
}

void kirkPatrick(List<Custom_Point> points) {

  print(Ordered.length);
  points.sort((a, b) => a.x.compareTo(b.x));
  Custom_Point pumin = findPumin(points);
  Custom_Point pumax = findPumax(points);
  Custom_Point plmin = findPlmin(points);
  Custom_Point plmax = findPlmax(points);
  print(pumin);
  print(pumax);
  print(plmin);
  print(plmax);
  Line medianLine = findMedianLine(points);

  List<Custom_Point> upperHullP = upperHull(pumin, pumax, points);

  List<Custom_Point> lowerHullP = lowerHull(plmin, plmax, points);
  /// Prepare the final convex hull by combining the upper and lower hulls and removing duplicates.
  upperHullP.sort((a, b) => a.x.compareTo(b.x));
  lowerHullP.sort((a, b) => b.x.compareTo(a.x));
  lowerHullP = lowerHullP.toSet().toList();
  upperHullP = upperHullP.toSet().toList();
  lowerHullP = removeConsecutiveDuplicates(lowerHullP);
  upperHullP = removeConsecutiveDuplicates(upperHullP);
  List <Custom_Point> convexHull = [];
  convexHull.addAll(upperHullP);
  convexHull.addAll(lowerHullP);
  convexHull.add(upperHullP[0]);
  // for(int i = 0; i < convexHull.length; i++) {
  //   fin.convexHull.add(convexHull[i]);
  // }
  
  // print("Convex Hull Points: ");
  // print(convexHull);
  cvhull= convexHull;
  print(Ordered.length);

}
