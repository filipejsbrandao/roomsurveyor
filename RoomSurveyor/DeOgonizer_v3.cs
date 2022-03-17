//
// DeOgonizer_v3.cs
//
// Author:
//       Filipe Jorge da Silva Brandao
//
// Copyright (c) 2020 ©2020 Filipe Brandao
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using RoomSurveyor;

namespace RoomSurveyor
{
    public class DeOgonizer_v3 : GH_Component
    {
        public DeOgonizer_v3()
          : base("DeOgonizer_v3", "DeOgon3",
            "This component transforms a ogon into a irregular non-convex polygon by recursively and randomly replacing square corners by non-orthogonal ones. Each corner is replaced by a random point inside a circle.",
            "RoomSurveyor", "Polygon Generators")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Ogon", "O", "A closed polyline that represents an orthogonal polygon", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Number", "N", "The number of corners to be replaced", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Seed", "S", "An integer to be used as a seed for the random generator. If nothing is provided the system clock is used. Provide a seed if you wish to debug your algorithm", GH_ParamAccess.item);
            pManager.AddNumberParameter("Ratio", "R", "A number between 0 and 1 that is used to scale the circle over which random points may be generated", GH_ParamAccess.item);
            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "The random polygon", GH_ParamAccess.item);
            pManager.AddPointParameter("Points", "Pt", "The new random corners of the polygon", GH_ParamAccess.list);
            pManager.AddCircleParameter("Circle", "C", "The circle over which the random points are generated", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Index", "I", "An ordered list of the indices of the changed corners", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double radius = 0.0;
            int n = 0;
            int s = 0;
            Polyline ogon = new Polyline();
            Curve curve = ogon.ToNurbsCurve();

            if (!DA.GetData(0, ref curve)) return;
            if (!DA.GetData(1, ref n)) return;
            bool user = DA.GetData(2, ref s);
            if (!DA.GetData(3, ref radius)) return;

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
            foreach (Point3d p in ogon)
            {
                //Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
                if (Math.Abs(first.Z - p.Z) > double.Epsilon)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The curve is not parallel to the XY plane. Consider providing the plane on which the curve lies or rotate it.");
                    break;
                }
            }

            if (n <= 0)
            { return; }
            if (radius > 1.0)
                radius = 1.0;
            if (radius <= 0.0)
                radius = 0.01;

            n = (n - 1) % ogon.SegmentCount;
            int count = ogon.SegmentCount;
            List<int> corners = new List<int>();

            Random r = new Random();
            Random rUser = new Random(s);

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

            List<Circle> circles = new List<Circle>();
            List<Point3d> pts = new List<Point3d>();
            List<Plane> pl = new List<Plane>();

            for (int i = 0; i < randomCorner.Count; i++)
            {
                int current = randomCorner[i];
                int next = (current + 1) % count;
                int prev = (current == 0) ? count - 1 : current - 1;
                Vector3d xAxis = new Vector3d(ogon[next] - ogon[prev]);
                Vector3d other = new Vector3d(ogon[current] - ogon[prev]);
                Plane plano = new Plane(ogon[prev], xAxis, other);
                plano.ClosestParameter(ogon[current], out double u, out double v);
                plano = new Plane(ogon[current], xAxis, other);
                pl.Add(plano);
                Circle circulo = new Circle(plano, v * radius);
                Vector3d raio = plano.YAxis * v * radius;
                var randRad = (!user) ? r.NextDouble() : rUser.NextDouble();
                var randY = (!user) ? r.NextDouble() : rUser.NextDouble();
                raio.Rotate(Math.PI * 2 * randRad, plano.ZAxis);
                Point3d newPoint = ogon[current] + raio * randY;
                circles.Add(circulo);
                pts.Add(newPoint);
                ogon.Insert(current, newPoint);
                ogon.RemoveAt(current + 1);
            }

            ogon.Add(ogon[0]);

            DA.SetData(0, ogon);
            DA.SetDataList(1, pts);
            DA.SetDataList(2, circles);
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
                return Properties.Resources.DeOgonizer_v3_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("2ca81dd1-ac76-472f-9e32-a4fde8aac913"); }
        }
    }
}
