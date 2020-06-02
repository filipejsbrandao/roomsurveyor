using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class OrientPolygon : GH_Component
    {
        public OrientPolygon()
          : base("Orient Polygon", "OrPoly",
            "Orients a polygon to counter-clockwise point order on a given plane, optionaly sets the start point of the polyline to the closest corner to a provided point",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polyline", "P", "A closed planar polyline. If the curve is not on the XY plane, please provide the plane the curve lies on", GH_ParamAccess.item);
            pManager.AddPlaneParameter("Plane", "Pl", "If the curve is not on the plane World XY provide the plane the curve lies on", GH_ParamAccess.item);
            pManager.AddPointParameter("Reference Point", "Pt", "A reference point to align the start point of the polyline to", GH_ParamAccess.item);
            pManager[1].Optional = true;
            pManager[2].Optional = true;
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
            Point3d refPt = Point3d.Origin;
            int shift;

            if (!DA.GetData(0, ref curve)) return;
            bool userPlane = DA.GetData(1, ref inPlane);
            bool userPt = DA.GetData(2, ref refPt);

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
                return;
            }

            //To ensure this always works we must move the polygon to the plane XY and then return it to it original position
            curve.TryGetPlane(out Plane curvePlane);

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

            if (!userPlane)
            {

                if (curvePlane.Normal.IsParallelTo(inPlane.Normal, 0.002) == 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The curve is not parallel to the XY plane within tolerance. Consider providing the plane on which the curve lies or rotate it.");
                    return;
                }
            }
            else
            {
                if (curvePlane.Normal.IsParallelTo(inPlane.Normal, 0.002) == 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon is not parallel to the provided plane");
                    return;
                }
            }

            Point3d startPt = Point3d.Origin;

            if (!userPlane)
            {
                poly = OrientPoly(poly, Plane.WorldXY);
                startPt = (userPt) ? poly[poly.ClosestIndex(refPt)] : poly[0];
                shift = (userPt) ? poly.ClosestIndex(refPt) : 0;

                if (startPt != poly[0])
                {
                    List<Point3d> points = new List<Point3d>();
                    points.AddRange(poly);
                    points.RemoveAt(points.Count - 1);
                    poly.Clear();
                    poly.AddRange(ShiftLeft(points, shift));
                    poly.Add(startPt);
                }
            }
            else
            {
                poly = OrientPoly(poly, inPlane);
                startPt = (userPt) ? poly[poly.ClosestIndex(refPt)] : poly[0];
                shift = (userPt) ? poly.ClosestIndex(refPt) : 0;

                if (startPt != poly[0])
                {
                    List<Point3d> points = new List<Point3d>();
                    points.AddRange(poly);
                    points.RemoveAt(points.Count - 1);
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
            if (orientation)
            {
                if (a.ClosedCurveOrientation(plane) == CurveOrientation.Clockwise)
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