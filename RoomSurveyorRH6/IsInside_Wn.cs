using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class IsInside_Wn : GH_Component
    {
        public IsInside_Wn()
          : base("Is Inside (Wn)", "IsIn",
            "Checks if some points are inside a polygon using the winding number counter algorithm. 0 = outside, 1 = the polygonal chain wraps around the point once, n = the polygon wraps n times around the point",
            "RoomSurveyor", "Ises")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "Pts", "A point to test inclusion", GH_ParamAccess.list);
            pManager.AddCurveParameter("Polygon", "P", "The Polygon to test the inclusion", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Winding Number", "Wn", "The number of times the polygon wraps around the point", GH_ParamAccess.list);
            pManager.AddPointParameter("Points Inside", "Pts'", "Points that are inside the curve projected onto the curve plane", GH_ParamAccess.list);
            pManager.AddPlaneParameter("Curve Plane", "Pl", "Plane the curve is on", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            List<Point3d> pts = new List<Point3d>();
            List<double> wnPoly = new List<double>(); // the number of times a polygon wraps arround a point
            List<Point3d> ptsInside = new List<Point3d>();

            if (!DA.GetDataList(0, pts)) return;
            if (!DA.GetData(1, ref curve)) return;

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The polygon must be a planar polyline");
                return;
            }

            //To ensure this always works we must move the polygon to the plane XY and then return it to it original position
            curve.TryGetPlane(out Plane uPlane);

            if (curve.TryGetPolyline(out poly))
            { }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "P input expects a polygon to be provided");
                return;
            }

            if (!poly.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Closed polylines must have at least 3 segments");
                return;
            }
            else if (!poly.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A closed polyline must be supplied");
                return;
            }
            foreach (Point3d pt in pts)
            {
                if (!pt.IsValid)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A valid point must be supplied");
                    return;
                }
            }

            Transform transform = Transform.PlaneToPlane(uPlane, Plane.WorldXY);
            poly.Transform(transform);
            int edges = poly.SegmentCount;
            Point3d[] pointArray = poly.ToArray();

            foreach (Point3d pt in pts)
            {
                Point3d tempPt = uPlane.ClosestPoint(pt);
                Point3d tempPt2 = new Point3d(tempPt);
                tempPt.Transform(transform);
                int wn = PnPoly(tempPt, pointArray, edges);
                wnPoly.Add(wn);
                if (wn != 0) ptsInside.Add(tempPt2);
            }

            DA.SetDataList(0, wnPoly);
            DA.SetDataList(1, ptsInside);
            DA.SetData(2, uPlane);

        }
        /// <summary>
        /// Winding Number test for a point in a polygon based on http://geomalgorithms.com/a03-_inclusion.html
        /// Copyright 2000 softSurfer, 2012 Dan Sunday
        /// This code may be freely used and modified for any purpose
        /// providing that this copyright notice is included with it.
        /// SoftSurfer makes no warranty for this code, and cannot be held
        /// liable for any real or imagined damage resulting from its use.
        /// Users of this code must verify correctness for their application.
        /// Ported to C# by Filipe Brandão in January 2019
        /// </summary>
        /// <returns> the winding number (=0 only when P is outside)</returns>
        /// <param name="Pt">Point.</param>
        /// <param name="poly">the Poly as an array of points where <paramref name="poly"/>[0] = poly[n]</param>
        /// <param name="n"> the number of edges in the polygon</param>
        public static int PnPoly(Point3d Pt, Point3d[] poly, int n)
        {
            int wn = 0; //the winding number counter

            for (int i = 0; i < n; i++)
            { //edge from V[i] to V[i+1]
                if (poly[i].Y <= Pt.Y)
                { // start Y <= Pt.Y
                    if (poly[i + 1].Y > Pt.Y) // an upward crossing
                        if (IsLeft(poly[i], poly[i + 1], Pt) > 0) // Pt is left of edge
                            ++wn; // have a valid up intersect
                }
                else
                { // start Y > Pt.Y (no test needed)
                    if (poly[i + 1].Y <= Pt.Y) // a downward crossing
                        if (IsLeft(poly[i], poly[i + 1], Pt) < 0) // Pt is right of edge
                            --wn; // have a valid down intersect
                }

            }
            return wn;
        }
        /// <summary>
        /// Determines if a third point is Left, On or Right of an infinite 2D line defined by two points
        /// </summary>
        /// <returns>Returns a positive number if it is Left, a negative number if it is Right and 0 if it is on the line</returns>
        /// <param name="P0">P0 - the first point</param>
        /// <param name="P1">P1 - the second point</param>
        /// <param name="P2">P2 - the point to test</param>

        public static double IsLeft(Point3d P0, Point3d P1, Point3d P2)
        {
            return ((P1.X - P0.X) * (P2.Y - P0.Y) - (P2.X - P0.X) * (P1.Y - P0.Y));
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.WnInside_Icon;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("39a742a8-3176-46c3-b3c8-5453051a5c85"); }
        }
    }
}
