class Custom_Point {
  final double x;
  final double y;

  Custom_Point(this.x, this.y);

  @override
  int get hashCode => x.hashCode ^ y.hashCode;

  @override
  String toString() => 'Custom_Point($x, $y)';
}

class animatedPoints {
  List<List<Custom_Point>> upperBridgePoints;
  List<List<Custom_Point>> lowerBridgePoints;
  List<List<Custom_Point>> quadrilateral;
  List<List<Custom_Point>> removedPoints;
  List<Custom_Point> final_convex_hull;
  List<List<List<Custom_Point>>> Ordered;

  animatedPoints(
    this.upperBridgePoints,
    this.lowerBridgePoints,
    this.quadrilateral,
    this.removedPoints,
      this.final_convex_hull,
      this.Ordered,
  );

  @override
  String toString() {
    String upperBridgePointsStr = upperBridgePoints.map((e) => e.toString()).join(', ');
    String lowerBridgePointsStr = lowerBridgePoints.map((e) => e.toString()).join(', ');
    String quadrilateralStr = quadrilateral.map((e) => e.toString()).join(', ');
    String removedPointsStr = removedPoints.map((e) => e.toString()).join(', ');

    return 'Animated Points:\n'
        'Upper Bridge Points: $upperBridgePointsStr\n'
        'Lower Bridge Points: $lowerBridgePointsStr\n'
        'Quadrilateral: $quadrilateralStr\n'
        'Removed Points: $removedPointsStr\n';
  }
}

class animatedPointsjm {
  List<List<Custom_Point>> updates;
  List<List<Custom_Point>> removedP;

  animatedPointsjm(
    this.updates,
    this.removedP,
  );

//   @override
//   String toString() {
//     String updates = updates.map((e) => e.toString()).join(', ');
//     String removedP = removedP.map((e) => e.toString()).join(', ');

//     return 'Animated Points:\n'
//         'Updates: $upperBridgePointsStr\n'
//         'Removed Points: $lowerBridgePointsStr\n'
//   }
}



