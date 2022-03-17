//
// ShiftStartPoint.cs
//
// Author:
//       Filipe Jorge da Silva Brandao
//
// Copyright (c) 2019 ©2019 Filipe Brandao
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

namespace RoomSurveyor

{
    public class ShiftStartPoint : GH_Component
    {

        public ShiftStartPoint()
          : base("ShiftStartPoint", "SSPt",
            "Shift the start point of a polygon by a number of steps",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A closed polyline or a polygon", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Steps", "S", "The number of steps to move the start point", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "Shifted polygon", GH_ParamAccess.item);
            pManager.AddPointParameter("Start Point", "SPt", "Start Point", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            int shift = 0;
            bool isCCW;
            bool direction;

            if (!DA.GetData(0, ref curve)) return;
            if (!DA.GetData(1, ref shift)) return;

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
                return;
            }

            curve.TryGetPlane(out Plane curvePlane);

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

            isCCW = PolyAngles.IsCCW(poly, curvePlane);
            direction = (shift > 0) ? true : false;
            if (!isCCW) direction = !direction;

            List<Point3d> points = new List<Point3d>();
            points.AddRange(poly);
            points.RemoveAt(points.Count - 1);
            shift = Math.Abs(shift) % points.Count;
            poly.Clear();
            if (direction)
            {
                poly.AddRange(OrientPolygon.ShiftLeft(points, shift));
                poly.Add(poly[0]);
            }
            else
            {
                poly.AddRange(OrientPolygon.ShiftRight(points, shift));
                poly.Add(poly[0]);
            }

            DA.SetData(0, poly);
            DA.SetData(1, poly[0]);
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.ShiftPolySPt_Icon;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("0311cf64-c883-4f96-9411-f9367945f231"); }
        }
    }
}
