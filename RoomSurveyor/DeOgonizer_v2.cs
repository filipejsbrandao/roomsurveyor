using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class DeOgonizer_v2 : GH_Component
    {
        public override Grasshopper.Kernel.GH_Exposure Exposure
        {
            get { return GH_Exposure.hidden; }
        }

        public DeOgonizer_v2()
          : base("DeOgonizer_v2", "DeOgon2",
            "This component transforms a ogon into a irregular non-convex polygon by recursively and randomly replacing square corners by non-orthogonal ones. Each corner is replaced by a random point inside a rectangle.",
            "RoomSurveyor", "Polygon Generators")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Ogon", "O", "A closed polyline that represents an orthogonal polygon", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Number", "N", "The number of corners to be replaced", GH_ParamAccess.item);
            pManager.AddNumberParameter("Seed", "S", "A number to be used as a seed for the random generator. If nothing is provided the system clock is used. Provide a seed if you wish to debug your algorithm", GH_ParamAccess.item);
            //pManager.AddNumberParameter("Ratio", "R", "A number between 0 and 1 that is used to scale the circle over which random points may be generated", GH_ParamAccess.item);
            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "The random polygon", GH_ParamAccess.item);
            pManager.AddPointParameter("Points", "Pt", "The new random corners of the polygon", GH_ParamAccess.list);
            pManager.AddRectangleParameter("Rectangle", "R", "The rectangle over which the random points are generated", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Index", "I", "An ordered list of the indices of the changed corners", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //double radius = 0.0;
            int n = 0;
            int seed = 0;
            Polyline ogon = new Polyline();
            Curve curve = ogon.ToNurbsCurve();

            if (!DA.GetData(0, ref curve)) return;
            if (!DA.GetData(1, ref n)) return;
            bool user = DA.GetData(2, ref seed);
            //if (!DA.GetData(3, ref radius)) return;

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
            if (!IsOrtho.IsPolyOrtho(ogon))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "BE WARNED: The provided polyline is not an orthogonal polygon, one or more of the internal angles are different from 90 or 270 degrees");
            }

            Point3d first = ogon[0];
            foreach (Point3d pt in ogon)
            {
                //Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
                if (Math.Abs(first.Z - pt.Z) > double.Epsilon)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The curve is not parallel to the XY plane. Consider providing the plane on which the curve lies or rotate it.");
                    break;
                }
            }

            if (n <= 0)
            { return; }

            n = (n - 1) % ogon.SegmentCount;
            int count = ogon.SegmentCount;
            List<int> corners = new List<int>();

            Random r = new Random();
            Random rUser = new Random((int)seed);

            List<int> randomCorner = new List<int>();

            ogon.RemoveAt(count);

            for (int i = 0; i < count; i++)
            {
                corners.Add(i);
            }

            for (int j = 0; j <= n; j++)
            {
                int rand = (!user) ? r.Next(corners.Count) : rUser.Next(corners.Count);
                randomCorner.Add(corners[rand]);
                corners.RemoveAt(rand);
            }

            List<Rectangle3d> p = new List<Rectangle3d>();
            List<Point3d> pts = new List<Point3d>();

            for (int i = 0; i < randomCorner.Count; i++)
            {
                int current = randomCorner[i];
                int next = (current + 1) % count;
                int prev = (current == 0) ? count - 1 : current - 1;
                Vector3d xAxis = new Vector3d(ogon[next] - ogon[prev]);
                Vector3d other = new Vector3d(ogon[current] - ogon[prev]);
                //Vector3d tolerance = xAxis * 0.1;
                Plane plano = new Plane(ogon[prev], xAxis, other);
                plano.ClosestParameter(ogon[current], out double u, out double v);
                plano.ClosestParameter(ogon[next], out double s, out double t);
                double rBase = (u < s - u) ? u : s - u;
                plano.Translate(plano.YAxis * v * 0.1 + plano.XAxis * (u - rBase));
                var randX = (!user) ? r.NextDouble() : rUser.NextDouble();
                var randY = (!user) ? r.NextDouble() : rUser.NextDouble();
                Point3d newPoint = plano.PointAt(randX * rBase * 2, randY * v * 1.5);
                Rectangle3d area = new Rectangle3d(plano, rBase * 2, v * 1.5);
                p.Add(area);
                pts.Add(newPoint);
                ogon.Insert(current, newPoint);
                ogon.RemoveAt(current + 1);

            }

            ogon.Add(ogon[0]);

            DA.SetData(0, ogon);
            DA.SetDataList(1, pts);
            DA.SetDataList(2, p);
            DA.SetDataList(3, randomCorner);
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.DeOgonizer_v2_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("471f6205-bf3a-4fdf-a01c-d9d55bf10d25"); }
        }
    }
}
