﻿//
// RoomTurtle.cs
//
// Author:
//       Filipe Jorge da Silva Brandao
//
// Copyright (c) 2021 ©2021 Filipe Brandao
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
    public class RoomTurtle : GH_Component
    {
        public RoomTurtle()
          : base("Room Turtle", "RT",
            "Creates an orthogonal polygonal chain based on the lengths of the sides and type of corner. True for PI/2 or 90 degree corner and false for 3/2 PI or 270 degree corner",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Side Length", "L", "The length of each side of the polygon in CCW order", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Angle", "A", "The type of angle at each corner. True for PI/2 or 90 degree corner and false for 3/2 PI or 270 degree corner", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Polyline", "P", "The polygonal chain is closed if the combination of suplied dimensions and angle terminate at the same point", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<double> length = new List<double>();
            List<bool> turn = new List<bool>();

            if ((length != null) && !DA.GetDataList(0, length)) return;
            if ((turn != null) && !DA.GetDataList(1, turn)) return;

            if (length.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A list of containing the lenght of each side must suplied");
                return;
            }
            if (turn.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A list of containing the a boolean representing the turn at each corner must suplied");
                return;
            }
            if (turn.Count < length.Count - 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The list of angles must have at least the one item less than the list of side lengths");
                return;
            }
            if (NumberCheck(length))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "All side lenghts must be positive real numbers");
                return;
            }

            int count = length.Count;
            Point3d currPt = Point3d.Origin;
            Vector3d nextV = new Vector3d(0, -1, 0);
            Polyline p = new Polyline();
            p.Add(currPt);

            for (int i = 0; i < count; i++)
            {
                currPt += nextV * length[i];

                p.Add(currPt);
                if (i == count - 1) { break; }
                if (turn[i])
                {
                    nextV.Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                }
                else
                {
                    nextV.Transform(Transform.Rotation(-Math.PI / 2, Point3d.Origin));
                }

            }
            DA.SetData(0, p);
        }

        

        /// <summary>
        /// Check if all numbers are positive and real.
        /// </summary>
        /// <returns><c>true</c>, if any number is not positive or real, <c>false</c> otherwise.</returns>
        /// <param name="a">The double numbers to check.</param>
        public static bool NumberCheck(List<double> a)
        {
            bool check = false;

            foreach (double n in a)
            {
                if (n <= 0)
                {
                    check = true;
                    break;
                }
                if (Double.IsNaN(n))
                {
                    check = true;
                    break;
                }
            }
            return check;
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.RoomTurtle_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("bbafe151-2002-4468-9ea1-66748094b9a5"); }
        }
    }
}
