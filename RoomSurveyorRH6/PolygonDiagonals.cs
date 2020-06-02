using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class PolygonDiagonals : GH_Component
    {
        public PolygonDiagonals()
          : base("Polygon Diagonals", "PolyDiag",
            "PolygonDiagonals description",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A closed planar polyline", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Internal Diagonals", "ID", "The internal diagonals of the polygon", GH_ParamAccess.tree);
            pManager.AddLineParameter("External Diagonals", "ED", "The diagonals that are outside of the polygon", GH_ParamAccess.tree);
            pManager.AddLineParameter("Intersecting Diagonals", "SID", "The diagonals of the polygon that intersect it", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            DataTree<Line> IntDiag = new DataTree<Line>();
            DataTree<Line> ExtDiag = new DataTree<Line>();
            DataTree<Line> SIDiag = new DataTree<Line>();

            if (!DA.GetData(0, ref curve)) return;

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
                return;
            }

            var events = Rhino.Geometry.Intersect.Intersection.CurveSelf(curve, 0.001);
            if (events.Count != 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The provided polygon is either self-intersecting or is very thin. Tolerance for self-intersections is 0.001");
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

            Transform transform = Transform.PlaneToPlane(curvePlane, Plane.WorldXY);
            poly.Transform(transform);

            transform = Transform.PlaneToPlane(Plane.WorldXY, curvePlane);

            int n = poly.Count - 1;

            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    if (i == 0)
                    {
                        if (j != i && j != i + 1 && j != poly.Count - 2)
                        {
                            GH_Path current = new GH_Path(i);
                            Line diagonal = new Line(poly[i], poly[j]);

                            int type = DiagonalType(diagonal, poly);
                            diagonal.Transform(transform);
                            switch (type)
                            {
                                case 1:
                                    IntDiag.Add(diagonal, current);
                                    break;
                                case -1:
                                    ExtDiag.Add(diagonal, current);
                                    break;
                                default:
                                    SIDiag.Add(diagonal, current);
                                    break;
                            }
                        }
                    }
                    else
                    {
                        if (j != i && j != i + 1)
                        {
                            GH_Path current = new GH_Path(i);
                            Line diagonal = new Line(poly[i], poly[j]);

                            int type = DiagonalType(diagonal, poly);
                            diagonal.Transform(transform);
                            switch (type)
                            {
                                case 1:
                                    IntDiag.Add(diagonal, current);
                                    break;
                                case -1:
                                    ExtDiag.Add(diagonal, current);
                                    break;
                                default:
                                    SIDiag.Add(diagonal, current);
                                    break;
                            }
                        }
                    }
                }
            }

            DA.SetDataTree(0, IntDiag);
            DA.SetDataTree(1, ExtDiag);
            DA.SetDataTree(2, SIDiag);
        }

        /// <summary>
        /// Returns the type of the diagonal.
        /// </summary>
        /// <returns><c>-1</c>an external diagonal.<c>1</c>an internal diagonal.<c>0</c> a self-intersecting diagonal.</returns>
        /// <param name="diagonal">Diagonal.</param>
        /// <param name="polyline">Polyline.</param>
        public static int DiagonalType(Line diagonal, Polyline polyline)
        {
            //a brute force method would be to test each edge of the polygon for intersection with the diagonal
            //if any intersection is found other than the start and end points the routine can be stopped where we can use Shamos-Hoey Algorithm
            //I think we could even simplify it by first testing one polygonal chain and then the other
            //we still need to test the diagonal for inclusion in the polygon
            LineCurve c1 = new LineCurve(diagonal);
            Curve c2 = polyline.ToNurbsCurve();

            const double intersection_tolerance = 0.001;
            const double overlap_tolerance = 0.0;
            var events = Rhino.Geometry.Intersect.Intersection.CurveCurve(c1, c2, intersection_tolerance, overlap_tolerance);
            int type = 0;

            if (events.Count == 2 && events[0].IsPoint && events[1].IsPoint)
            {
                //it has been garanteed elsewhere that the polyline is CCW oriented,
                //current is the Point at the origin of the diagonal which is also a point i of the polygon,
                //previous is the point [i-1] in the polygon and next is the point [i+1] in the polygon
                //we want to check if the diagonal formed by the current Point and the diagonalEnd bisects the angle
                Point3d current = diagonal.From;
                Point3d diagonalEnd = diagonal.To;
                int i = polyline.ClosestIndex(current);
                Point3d previous = Point3d.Unset;
                if (i == 0)
                {
                    previous = polyline[polyline.Count - 2];
                }
                else
                {
                    previous = polyline[i - 1];
                }
                Point3d next = polyline[i + 1];

                double first = RoomSurvey4.AngleAtCorner(current, previous, next); //must be Left so larger than 0
                double second = RoomSurvey4.AngleAtCorner(current, diagonalEnd, next); //must be Left so larger than 0

                if (first > second)
                    type = 1;
                else
                    type = -1;
            }

            return type;
        }
        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.PolyDiagonals2_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("9c9b769b-a945-4b44-9955-e82cf5e75393"); }
        }
    }
}
