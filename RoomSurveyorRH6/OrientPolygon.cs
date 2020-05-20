using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyorRH6
{
    public class OrientPolygon : GH_Component
    {
        public OrientPolygon()
          : base("Orient Polygon", "OrPoly",
            "A component that checks if an input polyline is a simple polygon, makes the polyline couterclockwise and sets the highest Y with lowest X as the start point",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polyline", "P", "A closed planar polyline. If the curve is not on the XY plane, please provide the plane the curve lies on", GH_ParamAccess.item);
            pManager.AddPlaneParameter("Plane", "O", "If the curve is not on the plane World XY provide the plane the curve lies on", GH_ParamAccess.item);
            pManager[1].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Oriented Polygon", "OP", "a counterclockwise oriented polygon", GH_ParamAccess.item);
            pManager.AddPointParameter("Start Point", "SPt", "the start point of the polyline", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            Plane inPlane = Plane.WorldXY;
            int shift;

            if (!DA.GetData(0, ref curve)) return;
            bool user = DA.GetData(1, ref inPlane);

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
                return;
            }

            if (curve.TryGetPolyline(out poly)) { }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "P input expects a polyline to be provided");
                return;
            }

            if (poly.IsValid == false)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Closed polylines must have at least 2 segments");
                return;
            }

            if (poly.IsClosed == false)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A closed polyline must be supplied");
                return;
            }

            Point3d first = poly[0];
            foreach (Point3d p in poly)
            {
                //Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
                if(Math.Abs(first.Z - p.Z) > double.Epsilon)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The curve is not parallel to the XY plane. Consider providing the plane on which the curve lies or rotate it.");
                    break;
                }
            }

            _ = Point3d.Origin;
            Point3d startPt;

            if (!user)
            {
                poly = OrientPoly(poly, Plane.WorldXY);
                startPt = FindStartPt(poly);

                if (startPt != poly[0])
                {
                    List<Point3d> points = new List<Point3d>();
                    points.AddRange(poly);
                    points.RemoveAt(points.Count - 1);
                    shift = StartPtIndex(poly);
                    poly.Clear();
                    poly.AddRange(ShiftLeft(points, shift));
                    poly.Add(startPt);
                }
            }
            else
            {
                poly = OrientPoly(poly, inPlane);
                startPt = FindStartPt(poly, inPlane);

                if (startPt != poly[0])
                {
                    List<Point3d> points = new List<Point3d>();
                    points.AddRange(poly);
                    points.RemoveAt(points.Count - 1);
                    shift = StartPtIndex(poly, inPlane);
                    poly.Clear();
                    poly.AddRange(ShiftLeft(points, shift));
                    poly.Add(startPt);
                }
            }

            DA.SetData(0, poly);
            DA.SetData(1, startPt);
        }

        /// <summary>
        /// Orient a closed polyline on a plane to a given orientation.
        /// </summary>
        /// <param name="poly">The polyline to orient</param>
        /// <param name="plane">The reference plane</param>
        /// <param name="orientation">if true it will orient the polyline in counter-clockwise orientation, false will return the polyline clockwise oriented</param>
        /// <returns></returns>
        public static Polyline OrientPoly(Polyline poly, Plane plane, bool orientation)
        {
            Curve a = poly.ToNurbsCurve();
            if (orientation) {
                if(a.ClosedCurveOrientation(plane) == CurveOrientation.Clockwise)
                poly.Reverse();
            }
            else
            {
                if (a.ClosedCurveOrientation(plane) == CurveOrientation.CounterClockwise)
                    poly.Reverse();
            }
            return poly;
        }

        /// <summary>
        /// Orient a closed polyline on a plane to a counter-clockwise orientation.
        /// </summary>
        /// <param name="poly">The polyline to orient</param>
        /// <param name="plane">The reference plane</param>
        /// <returns></returns>
        public static Polyline OrientPoly(Polyline poly, Plane plane)
        {
            Curve a = poly.ToNurbsCurve();
            if (a.ClosedCurveOrientation(plane) == CurveOrientation.Clockwise)
                poly.Reverse();

            return poly;
        }

        /// <summary>
        /// Orient a closed polyline on the XY plane to counter-clockwise orientation.
        /// </summary>
        /// <param name="poly"></param>
        /// <returns></returns>
        public static Polyline OrientPoly(Polyline poly)
        {
            Vector3d unitZ = new Vector3d(0, 0, 1);
            Curve a = poly.ToNurbsCurve();
            if (a.ClosedCurveOrientation(unitZ) == CurveOrientation.Clockwise)
                poly.Reverse();

            return poly;
        }

        /// <summary>
        /// Find the index of corner of the polygon with largest Y and smallest X coordinates
        /// </summary>
        /// <returns>The point index.</returns>
        /// <param name="poly">Poly.</param>
        public static int StartPtIndex(Polyline poly)
        {

            int index = 0;
            Point3d pt = poly[index];

            for (int i = 1; i < poly.Count; i++)
            {
                if (poly[i].Y > pt.Y)
                {
                    index = i;
                    pt = poly[i];
                }
                else if (Math.Abs(poly[i].Y - pt.Y) < double.Epsilon && poly[i].X < pt.X)
                {
                    index = i;
                    pt = poly[i];
                }
            }
            return index;
        }

        /// <summary>
        /// Find the index of corner of the polygon with largest Y and smallest X coordinates
        /// </summary>
        /// <returns>The point index.</returns>
        /// <param name="poly">Poly.</param>
        public static int StartPtIndex(Polyline poly, Plane plane)
        {

            int index = 0;
            plane.RemapToPlaneSpace(poly[index], out Point3d pt);

            for (int i = 1; i < poly.Count; i++)
            {
                plane.RemapToPlaneSpace(poly[i], out Point3d planePt);
                if (planePt.Y > pt.Y)
                {
                    index = i;
                    pt = planePt;
                }
                else if (Math.Abs(planePt.Y - pt.Y) < double.Epsilon && planePt.X < pt.X)
                {
                    index = i;
                    pt = planePt;
                }
            }
            return index;
        }

        /// <summary>
        /// Find the corner of the polygon with largest Y and smallest X coordinates
        /// </summary>
        /// <returns>The start point.</returns>
        /// <param name="poly">Poly.</param>
        public static Point3d FindStartPt(Polyline poly)
        {

            Point3d pt = poly[0];

            foreach (Point3d p in poly)
                if (p.Y > pt.Y)
                    pt = p;
                else if (Math.Abs(p.Y - pt.Y) < Double.Epsilon && p.X < pt.X)
                    pt = p;

            return pt;
        }

        /// <summary>
        /// Find the corner of the polygon with largest Y and smallest X coordinates
        /// </summary>
        /// <returns>The start point.</returns>
        /// <param name="poly">Poly.</param>
        public static Point3d FindStartPt(Polyline poly, Plane plane)
        {
            plane.RemapToPlaneSpace(poly[0], out Point3d pt);
            int i = 0;
            int startIndex = 0;

            foreach (Point3d p in poly)
            {
                plane.RemapToPlaneSpace(p, out Point3d planePt);
                if (planePt.Y > pt.Y)
                {
                    pt = planePt;
                    startIndex = i;
                }
                else if (Math.Abs(planePt.Y - pt.Y) < Double.Epsilon && planePt.X < pt.X) {
                    pt = planePt;
                    startIndex = i;
                }
                    
                i++;
            }

            return poly[startIndex];
        }

        /// <summary>
        /// Shifts a list of points to the left or downwards
        /// </summary>
        /// <returns>The  shifted list to the left.</returns>
        /// <param name="list">List of points.</param>
        /// <param name="shiftBy">The number to shift the list by.</param>
        public static List<Point3d> ShiftLeft(List<Point3d> list, int shiftBy)
        {
            if (list.Count <= shiftBy)
            {
                return list;
            }

            var result = list.GetRange(shiftBy, list.Count - shiftBy);
            result.AddRange(list.GetRange(0, shiftBy));

            return result;
        }

        /// <summary>
        /// Shifts a list of points to the right or upwards
        /// </summary>
        /// <returns>The shifted list.</returns>
        /// <param name="list">a List of points.</param>
        /// <param name="shiftBy">The number to shift the list by.</param>
        public static List<Point3d> ShiftRight(List<Point3d> list, int shiftBy)
        {
            if (list.Count <= shiftBy)
            {
                return list;
            }

            var result = list.GetRange(list.Count - shiftBy, shiftBy);
            result.AddRange(list.GetRange(0, list.Count - shiftBy));
            return result;
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.OrientPoly_Icon;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("e95bd252-c2e3-41d4-b353-b9f7787c37b7"); }
        }
    }
}