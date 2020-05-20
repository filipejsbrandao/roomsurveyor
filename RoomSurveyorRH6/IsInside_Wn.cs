using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyorRH6
{
    public class IsInside_Wn : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public IsInside_Wn()
          : base("Is Inside (Wn)", "IsIn",
            "Checks if a point is inside a polygon using the winding number counter algorithm. If the polygon is oriented in the Counter Clockwise direction the algorithm returns positive numbers, otherwise it returns negative numbers for inclusion. 0 = outside, 1 = the polygonal chain wraps around the point once, n = the polygon wraps n times around the point",
            "RoomSurveyor", "Ises")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Point", "Pt", "A point to test inclusion", GH_ParamAccess.item);
            pManager.AddCurveParameter("Polygon", "P", "The Polygon to test the inclusion", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Winding Number", "WN", "The number of times the polygon wraps around the point", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            Point3d pt = new Point3d();
            double wnPoly; // the number of times a polygon wraps arround a point

            if (!DA.GetData(0, ref pt)) return;
            if (!DA.GetData(1, ref curve)) return;

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The polygon must be a planar polyline");
                return;
            }

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
            else if (!pt.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A valid point must be supplied");
                return;
            }

            int edges = poly.SegmentCount;
            Point3d[] pointArray = poly.ToArray();

            wnPoly = PnPoly(pt, pointArray, edges);

            DA.SetData(0, wnPoly);

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
