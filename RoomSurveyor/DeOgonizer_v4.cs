using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class DeOgonizer_v4 : GH_Component
    {

        public DeOgonizer_v4()
          : base("DeOgonizer_v4", "DeOgon4",
            "This algorithm randomly transforms an ogon into an non-orthogonal polygon. It does it by randomly moving a corner along the infinite line of the previous edge line segment, back or forth by half the length of the segment. It implements a few methods to prevent transformations that may generate self-intersecting polygons.",
            "RoomSurveyor", "Polygon Generators")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Ogon", "O", "A closed polyline that represents an orthogonal polygon", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Number", "N", "The number of corners to be replaced", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Seed", "S", "An integer to be used as a seed for the random generator. If nothing is provided the system clock is used. Provide a seed if you wish to debug your algorithm", GH_ParamAccess.item);
            pManager.AddNumberParameter("Ratio", "R", "A number between 0 and 1 that is used to scale the circle over which random points may be generated", GH_ParamAccess.item, 0.5);
            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "The random polygon", GH_ParamAccess.item);
            pManager.AddPointParameter("Points", "Pt", "The new random corners of the polygon", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Index", "I", "An ordered list of the indices of the changed corners", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double ratio = 0.0;
            int n = 0;
            int s = 0;
            Polyline ogon = new Polyline();
            Curve curve = ogon.ToNurbsCurve();

            if (!DA.GetData(0, ref curve)) return;
            if (!DA.GetData(1, ref n)) return;
            bool user = DA.GetData(2, ref s);
            if (!DA.GetData(3, ref ratio)) return;

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

            if (curve.TryGetPolyline(out ogon))
            { }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "P input expects a polygon to be provided");
                return;
            }

            if (!ogon.IsValid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Closed polylines must have at least 3 segments");
                return;
            }
            if (!ogon.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A closed polyline must be supplied");
                return;
            }

            if (n <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "N must be larger than 0");
                return;
            }

            if (!IsOrtho.IsPolyOrtho(ogon))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "BE WARNED: The provided polyline is not an orthogonal polygon, one or more of the internal angles are different from 90 or 270 degrees");
            }

            Point3d first = ogon[0];
            foreach (Point3d p in ogon)
            {
                //Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
                if (Math.Abs(first.Z - p.Z) > double.Epsilon)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The curve is not parallel to the XY plane. Consider providing the plane on which the curve lies or rotate it.");
                    break;
                }
            }

            if (ratio > 1.0)
                ratio = 1.0;
            if (ratio <= 0.0)
                ratio = 0.01;



            //Here we make sure that the number of corners changed is smaller or equal to the number of corners of the polygon
            n = (n - 1) % ogon.SegmentCount;
            int count = ogon.SegmentCount;
            List<int> corners = new List<int>();

            Random r = new Random();
            Random rUser = new Random(s);
            List<int> randomCorner = new List<int>();

            //here we create a list of corners indices
            for (int i = 0; i < count; i++)
            {
                corners.Add(i);
            }

            //here we randomly select one of said corners and build a list with it
            for (int j = 0; j <= n; j++)
            {
                int rand = (!user) ? r.Next(corners.Count) : rUser.Next(corners.Count);
                randomCorner.Add(corners[rand]);
                corners.RemoveAt(rand);
            }

            List<Point3d> pts = new List<Point3d>();

            //This does all the work
            for (int i = 0; i < randomCorner.Count; i++)
            {
                int current = randomCorner[i];
                int next = (current + 1) % count;
                int prev = (current == 0) ? count - 1 : current - 1;
                var randL = (!user) ? r.NextDouble() : rUser.NextDouble();
                double halfL = (ogon.SegmentAt(prev).Length / 2.0) * ratio;
                Vector3d move = ogon[current] - ogon[prev];
                move.Unitize();
                int turn = (!user) ? r.Next(2) : rUser.Next(2);
                move = (turn < 1) ? move * halfL * randL : -move * halfL * randL;
                Point3d newPoint = ogon[current] + move;
                if (SelfIntersect(ogon[current], ogon[next], newPoint, ogon, out double ccx_distance))
                {
                    int other = ClosestCorner(ogon, ogon[current], ogon[next], newPoint, out double distance);
                    if (ccx_distance < distance || other < 0)
                        distance = ccx_distance;
                    distance /= 3.0;
                    Vector3d move2 = move;
                    move2.Unitize();
                    newPoint = newPoint + (-1 * move) + move2 * distance;
                }
                pts.Add(newPoint);
                ogon.Insert(current, newPoint);
                ogon.RemoveAt(current + 1);
                if (current == 0)
                {
                    ogon.RemoveAt(ogon.SegmentCount);
                    ogon.Add(ogon[0]);
                }
            }

            DA.SetData(0, ogon);
            DA.SetDataList(1, pts);
            DA.SetDataList(2, randomCorner);
        }

        /// <summary>
        /// Determines if a line segment between two points intersects a polygon. The first point must be a corner of the polygon.
        /// </summary>
        /// <param name="current">The current corner of the polygon</param>
        /// <param name="next">The next corner of the polygon</param>
        /// <param name="newPoint">The new corner of the polygon</param>
        /// <param name="polygon">a closed polyline</param>
        /// <param name="distance">returns the distance of the closest intersection with the polygon to the infinite line from current to next .</param>
        /// <returns>true if the line intersects the polygon</returns>
        public bool SelfIntersect(Point3d current, Point3d next, Point3d newPoint, Polyline polygon, out double distance)
        {
            bool intersect = false;
            Curve c1 = polygon.ToNurbsCurve();
            const double intersection_tolerance = 0.001;
            const double overlap_tolerance = 0.0;
            Curve c2 = new Polyline(new List<Point3d>() { current, newPoint, next }).ToNurbsCurve();
            Line edge = new Line(current, next);

            var events = Rhino.Geometry.Intersect.Intersection.CurveCurve(c1, c2, intersection_tolerance, overlap_tolerance);
            distance = edge.DistanceTo(newPoint, false);

            if (events.Count > 2)
            {
                intersect = true;

                for (int i = 0; i < events.Count; i++)
                {
                    var ccx_event = events[i];
                    if (ccx_event.PointA.EpsilonEquals(current, double.Epsilon))
                        continue;
                    if (edge.DistanceTo(ccx_event.PointA, false) < distance)
                        distance = edge.DistanceTo(ccx_event.PointA, false);
                }
            }

            return intersect;
        }

        /// <summary>
        /// This methods checks if any corner of the polygon is inside a given triangle and returns the index of the closest corner of the polygon
        /// </summary>
        /// <param name="polygon">A polyline representing the polygon</param>
        /// <param name="current">The current corner of the polygon</param>
        /// <param name="next">The next corner of the polygon</param>
        /// <param name="newPoint">The point to be added to the polygon</param>
        /// <returns>The index of the closest corner of the polygon to the line from current to next or -1 if no corner lies inside the triangle</returns>
        public int ClosestCorner(Polyline polygon, Point3d current, Point3d next, Point3d newPoint, out double distance)
        {
            int closest_corner = -1;
            Line ln = new Line(current, next);
            distance = ln.DistanceTo(newPoint, false);

            Point3d[] triangle = { current, next, newPoint, current };

            for (int i = 0; i < polygon.SegmentCount; i++)
            {
                if (Math.Abs(IsInside_Wn.PnPoly(polygon[i], triangle, 3)) > 0)
                {
                    if (distance > ln.DistanceTo(polygon[i], true))
                    {
                        closest_corner = i;
                        distance = ln.DistanceTo(polygon[i], true);
                    }
                }
            }
            return closest_corner;
        }
        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.DeOgonizer_v4_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("fec9f1a1-e293-4511-ba91-0bb707261b8a"); }
        }
    }
}
