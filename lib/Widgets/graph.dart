import 'package:flutter/widgets.dart';

import 'kps.dart';
import 'jm.dart';
import 'package:flutter/material.dart';
import 'package:syncfusion_flutter_charts/charts.dart';
import 'Custom_Point.dart';

class Graph extends StatefulWidget {
  const Graph({super.key});

  @override
  State<Graph> createState() => _GraphState();
}

class _GraphState extends State<Graph> {

  List<Custom_Point> data = [
    Custom_Point(30, 70),
    Custom_Point(53, 32),
    Custom_Point(74, 66),
    Custom_Point(40, 80),
    Custom_Point(17, 25),
    Custom_Point(80, 45),
    Custom_Point(10, 90),
    Custom_Point(60, 10),
    Custom_Point(90, 20),
    Custom_Point(45, 55),
    Custom_Point(85, 85),
    Custom_Point(25, 35),
    Custom_Point(15, 60),
    Custom_Point(70, 5),
    Custom_Point(5, 10),
    Custom_Point(50, 75),
    Custom_Point(95, 50),
    Custom_Point(20, 30),
    Custom_Point(35, 65),
    Custom_Point(75, 15),
    Custom_Point(55, 80),
    Custom_Point(65, 40),
    Custom_Point(45, 20),
    Custom_Point(10, 40),
    Custom_Point(80, 75),
    Custom_Point(30, 50),
    Custom_Point(55, 5),
    Custom_Point(85, 30),
    Custom_Point(25, 75),
    Custom_Point(60, 90),
    Custom_Point(40, 10),
    Custom_Point(90, 60),
    Custom_Point(70, 35),
    Custom_Point(50, 30),
    Custom_Point(15, 5),
    Custom_Point(20, 45),
    Custom_Point(75, 80),
    Custom_Point(5, 85),
    Custom_Point(35, 15),
    Custom_Point(95, 10),
    Custom_Point(10, 65),
    Custom_Point(45, 50),
    Custom_Point(60, 25),
    Custom_Point(25, 10),
    Custom_Point(50, 60),
    Custom_Point(80, 20),
    Custom_Point(30, 95),
    Custom_Point(55, 70),
    Custom_Point(85, 40),
    Custom_Point(70, 85),
    Custom_Point(40, 55),
    Custom_Point(15, 20),
    Custom_Point(5, 50),
    Custom_Point(45, 5),
    Custom_Point(90, 90),
    Custom_Point(20, 80),
    Custom_Point(75, 25),
    Custom_Point(35, 45),
    Custom_Point(95, 75),
    Custom_Point(10, 15),
    Custom_Point(65, 80),
    Custom_Point(25, 60),
    Custom_Point(60, 35),
    Custom_Point(50, 10),
    Custom_Point(80, 65),
    Custom_Point(30, 25),
    Custom_Point(55, 90),
    Custom_Point(85, 15),
    Custom_Point(70, 50),
    Custom_Point(40, 40),
    Custom_Point(15, 75),
    Custom_Point(5, 20),
    Custom_Point(45, 95),
    Custom_Point(90, 30),
    Custom_Point(20, 5),
    Custom_Point(75, 60),
    Custom_Point(35, 80),
    Custom_Point(95, 20),
    Custom_Point(10, 50),
    Custom_Point(65, 15),
    Custom_Point(25, 90),
    Custom_Point(60, 55),
    Custom_Point(50, 35),
    Custom_Point(85, 70),
    Custom_Point(70, 95),
    Custom_Point(40, 25),
    Custom_Point(15, 10),
    Custom_Point(5, 65),
    Custom_Point(45, 30),
    Custom_Point(90, 5),
    Custom_Point(20, 50),
    Custom_Point(75, 85),
    Custom_Point(35, 10),
    Custom_Point(95, 45),
    Custom_Point(10, 80),
    Custom_Point(65, 25),
    Custom_Point(25, 5),
    Custom_Point(60, 70),
    Custom_Point(50, 45),
    Custom_Point(85, 10),
    Custom_Point(70, 40),
    Custom_Point(40, 95),
    Custom_Point(15, 50),
    Custom_Point(5, 85),
    Custom_Point(45, 15),
    Custom_Point(90, 60),
    Custom_Point(20, 35),
    Custom_Point(75, 10),
    Custom_Point(35, 60),
    Custom_Point(95, 30)
  ];
  List<Custom_Point> final_convex_hull = [];

  List<List<Custom_Point>> upperBridgePoints = [];
  List<List<Custom_Point>> lowerBridgePoints = [];
  List<Custom_Point> quadrilateral = [];
  List<List<Custom_Point>> removedPoints = [];
  List<List<List<Custom_Point>>> Ordered = [];
  int splice = 0;
  Custom_Point p = Custom_Point(0, 0);
  Custom_Point q = Custom_Point(0, 0);
  List<Custom_Point> currline = [Custom_Point(0, 0), Custom_Point(0, 0)];

  List<List<Custom_Point>> updates = [];
  List<Custom_Point> hull = [];
  List<List<Custom_Point>> removedP = [];
  List<Custom_Point> temp = [];

  double pointx = 0;
  double pointy = 0;

  void addPointx(String point) {
    setState(() {
      pointx = double.parse(point);
    });
  }

  void addPointy(String point) {
    setState(() {
      pointy = double.parse(point);
    });
  }

  void addPoint() {
    setState(() {
      data.add(Custom_Point(pointx, pointy));
    });
  }

  void runAlgo() async {
    animatedPoints ap = kirkPatrick(data);
    print(ap);
    setState(() {
      p=q;
      currline=[];
      updates = [];
      hull = [];
      removedP = [];
      temp = [];
      Ordered = ap.Ordered;
      // final_convex_hull = ap.final_convex_hull;
    });
    int a = 10;
    for (int i = 0; i < Ordered.length; i++) {
      setState(() {
        upperBridgePoints.add(Ordered[i][0]);
      });

      await Future.delayed(Duration(milliseconds: a));
      setState(() {
        quadrilateral = Ordered[i][1];
      });
      if(Ordered[i][1].isEmpty) {
        continue;
      }

      await Future.delayed(const Duration(milliseconds: 1000));
      setState(() {
        removedPoints.add(Ordered[i][2]);
      });

      await Future.delayed(const Duration(milliseconds: 1000));
    }
    await Future.delayed(const Duration(milliseconds: 1000));
    setState(() {
  final_convex_hull = ap.final_convex_hull;
});
  }

  void runAlgojm() async {
    
    animatedPointsjm jm = jarvisMarch(data);

    setState(() {
   final_convex_hull = [];

 upperBridgePoints = [];
 lowerBridgePoints = [];
 quadrilateral = [];
 removedPoints = [];
Ordered = [];
      updates = jm.updates;
      removedP = jm.removedP;
      print(updates);
      hull.add(updates[0][0]);
    });

    for (int i = 0; i < updates.length; i++) {
      setState(() {
        if (i >= 3) {
          temp += removedP[i - 3];
        }
        p = updates[i][0];
        currline[0] = p;
      });

      for (int j = 1; j < updates[i].length; j++) {
        setState(() {
          q = updates[i][j];
          currline[1] = q;
        });
        await Future.delayed(const Duration(milliseconds: 300));
      }

      await Future.delayed(const Duration(milliseconds: 600));
      setState(() {
        hull.add(updates[i][updates[i].length - 1]);
      });
    }
    setState(() {
      hull.add(hull[0]);
      p = q;
    });
  }

  void clearPoints() {
    setState(() {
      data = [];
      final_convex_hull = [];
      upperBridgePoints = [];
      lowerBridgePoints = [];
      quadrilateral = [];
      removedPoints = [];
      Ordered = [];
      p=q;
      currline=[];
      updates = [];
      hull = [];
      removedP = [];
      temp = [];
    });
  }

  @override
  Widget build(BuildContext context) {
    return SingleChildScrollView(
      child:Padding(padding: EdgeInsets.all(20),
      child: Column(
        children: [
          SfCartesianChart(
            margin: EdgeInsets.all(20),
            primaryXAxis: const NumericAxis(
              
              minimum: -10,
              maximum: 100,
                title: AxisTitle(text: 'X'),
                crossesAt: 0,
                placeLabelsNearAxisLine: false),
            primaryYAxis: const NumericAxis(
              minimum: -10,
              maximum: 100,
                title: AxisTitle(text: 'Y'),
                crossesAt: 0,
                placeLabelsNearAxisLine: false),
            series: <CartesianSeries>[
              ScatterSeries<Custom_Point, double>(
                dataSource: data,
                xValueMapper: (Custom_Point data, _) => data.x as double?,
                yValueMapper: (Custom_Point data, _) => data.y,
              ),
              ...upperBridgePoints.map((List<Custom_Point> points) {
                return LineSeries<Custom_Point, double>(
                    animationDuration: 0,
                    
                    dataSource: points,
                    xValueMapper: (Custom_Point data, _) => data.x as double?,
                    yValueMapper: (Custom_Point data, _) => data.y,
                    color: Colors.green,
                    width: 2);
              }),
              LineSeries<Custom_Point, double>(
                animationDuration: 0,
               
                dataSource: quadrilateral,
                xValueMapper: (Custom_Point data, _) => data.x as double?,
                yValueMapper: (Custom_Point data, _) => data.y,
                color: Colors.red,
              ),
              ...removedPoints.map((List<Custom_Point> points) {
                return ScatterSeries<Custom_Point, double>(
                  animationDuration: 0,
                  
                  dataSource: points,
                  xValueMapper: (Custom_Point data, _) => data.x as double?,
                  yValueMapper: (Custom_Point data, _) => data.y,
                  color: Colors.red,
                  // pointColorMapper: (Custom_Point data, _) => Colors.blue,
                  // pointBorderColorMapper: (Custom_Point data, _) => Colors.blue,
                  // borderWidth: 2
                );
              }),

              LineSeries<Custom_Point, double>(
                  animationDuration: 0,
                  dataSource: final_convex_hull,
                  xValueMapper: (Custom_Point data, _) => data.x as double?,
                  yValueMapper: (Custom_Point data, _) => data.y,
                  color: Colors.blue,
                  width: 2),
              LineSeries<Custom_Point, double>(
                  animationDuration: 100,
                  dataSource: hull,
                  xValueMapper: (Custom_Point data, _) => data.x as double?,
                  yValueMapper: (Custom_Point data, _) => data.y,
                  color: Colors.green,
                  width: 2),
              LineSeries<Custom_Point, double>(
                  animationDuration: 0,
                  dataSource: [p, q],
                  xValueMapper: (Custom_Point data, _) => data.x as double?,
                  yValueMapper: (Custom_Point data, _) => data.y,
                  color: Colors.red,
                  width: 2),
              ScatterSeries<Custom_Point, double>(
                animationDuration: 0,
                dataSource: temp,
                xValueMapper: (Custom_Point data, _) => data.x as double?,
                yValueMapper: (Custom_Point data, _) => data.y,
                color: Colors.red,
              ),
            ],
          ),
          Container(
            padding: const EdgeInsets.fromLTRB(500, 50, 500, 10),
            child: TextField(
              decoration: const InputDecoration(
                border: OutlineInputBorder(),
                hintText: 'X',
              ),
              onChanged: (String value) {
                addPointx(value);
              },
            ),
          ),
          Container(
            padding: const EdgeInsets.fromLTRB(500, 10, 500, 20),
            child: TextField(
              decoration: const InputDecoration(
                border: OutlineInputBorder(),
                hintText: 'Y',
              ),
              onChanged: (String value) {
                addPointy(value);
              },
            ),
          ),
          const SizedBox(width: 20),
          Row(
            mainAxisAlignment: MainAxisAlignment.center,
            children: [
              ElevatedButton(
                onPressed: addPoint,
                child: const Text('Add Custom_Point'),
              ),
              ElevatedButton(
                onPressed: clearPoints,
                child: const Text('Clear Points'),
              ),
            ],
          ),
          Row(
            mainAxisAlignment: MainAxisAlignment.center,
            children: [
              ElevatedButton(
                onPressed: runAlgo,
                child: const Text('Run kps'),
              ),
              ElevatedButton(
                onPressed: () {
                  setState(() {
                   
                  });
                  runAlgojm();},
                child: const Text('Run jm'),
              ),
            ],
          ),
        ],
      )),
    );
  }
}
