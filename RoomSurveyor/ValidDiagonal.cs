//
// ValidDiagonal.cs
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
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class ValidDiagonal : GH_Component
    {
        public ValidDiagonal()
          : base("Internal Diagonals", "ID",
            "This component returns all internal diagonals of a polygon for each corner without duplicates as a tree. A tree of the indexes of the diagonal endpoints for each corner is also provided.",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A closed planar polygon", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Diagonals", "D", "A tree with the internal diagonals of the polygon, duplicates are removed", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Indices", "I", "Returns the indices of the other endpoints of each of the internal diagonals for each corner", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {

            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();

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


            //Isto foi comentado para testar um problema com as trees e listas
            //Ou seja o codigo funciona apenas se for dado um elemento em todos os outros casos junta todos os elementos
            //numa unica árvore (tree) em que os indices são os cantos do poligono
            
            Transform transform = Transform.PlaneToPlane(curvePlane, Plane.WorldXY);
            poly.Transform(transform);

            transform = Transform.PlaneToPlane(Plane.WorldXY, curvePlane);
            double[,] diagMatrix = Triangulation.UpperMatrix(poly);
            //var diagonals = new GH_Structure<GH_Line>();
            //var indices = new GH_Structure<GH_Integer>();
            //diagonals.AppendRange(DiagonalTree1(diagMatrix, poly));
            //indices.AppendRange(IndexTree1(diagMatrix));
            var diagonals = DiagonalTree1(diagMatrix, poly);
            var indices = IndexTree1(diagMatrix);

            //This needs to be checked if it works (if it translates the diagonals back to the original polygon positions)
            for (int pi = 0; pi < diagonals.PathCount; pi++)
            {
                var path = diagonals.Paths[pi];
                var list = diagonals.Branches[pi];
                for (int i = 0; i < list.Count; i++)
                    if (list[i] != null)
                        if (list[i].Value != null)
                        {
                            Line l = list[i].Value;
                            l.Transform(transform);
                            list[i].Value = l;

                        }

                //collection.Add(list[i].Value);
            }
            /*
            foreach (var path in diagonals.Paths)
            {
                //List<Line> linesAtPath = diagonals.Branchpath);
                //List<GH_Line> linesAtPath = diagonals.get_Branch(path);
                for (int i = 0; i < linesAtPath.Count; i++)
                {
                    Line l = linesAtPath[i];
                    l.Transform(transform);
                    linesAtPath[i] = l;
                }
            }*/

            DA.SetDataTree(0, diagonals);
            DA.SetDataTree(1, indices);

        }

        /// <summary>
        /// Creates a tree of lines of the valid diagonals of the polygon.
        /// </summary>
        /// <param name="diagMatrix">the matrix containing the diagonals of the polygon.</param>
        /// <returns>Returns a tree of lines containing the diagonals of the polygon without duplicates.</returns>
        public GH_Structure<GH_Line> DiagonalTree1(double[,] diagMatrix, Polyline poly)
        {
            int matrixSize = diagMatrix.GetLength(0);
            int kCount = matrixSize * (matrixSize + 1) / 2;
            GH_Structure<GH_Line> diagonals = new GH_Structure<GH_Line>();
            for (int k = 0; k < kCount; k++)
            {
                //to access the UpperMatrix
                int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                int j = k + i * (i + 1) / 2 - matrixSize * i;
                if (diagMatrix[i, j] > 0)
                {
                    GH_Path current = new GH_Path(i);
                    diagonals.Append(new GH_Line(new Line(poly[i], poly[j])), current);                    
                }
            }
            return diagonals;
        }

        /// <summary>
        /// Creates a tree of lines of the valid diagonals of the polygon.
        /// </summary>
        /// <param name="diagMatrix">the matrix containing the diagonals of the polygon.</param>
        /// <returns>Returns a tree of lines containing the diagonals of the polygon without duplicates.</returns>
        public DataTree<Line> DiagonalTree(double[,] diagMatrix, Polyline poly)
        {
            int matrixSize = diagMatrix.GetLength(0);
            int kCount = matrixSize * (matrixSize + 1) / 2;
            DataTree<Line> diagonals = new DataTree<Line>();
            for (int k = 0; k < kCount; k++)
            {
                //to access the UpperMatrix
                int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                int j = k + i * (i + 1) / 2 - matrixSize * i;
                if (diagMatrix[i, j] > 0)
                {
                    GH_Path current = new GH_Path(i);
                    diagonals.Add(new Line(poly[i], poly[j]), current);
                }
            }

            return diagonals;
        }

        /// <summary>
        /// Creates a tree of the indices of the opposite corners of the valid diagonals of the polygon that connect each corner.
        /// </summary>
        /// <param name="diagMatrix">the matrix containing the diagonals of the polygon</param>
        /// <returns>Returns a tree of the indices of the opposite corner of each diagonal that connect to each corner.</returns>
        public GH_Structure<GH_Integer> IndexTree1(double[,] diagMatrix)
        {
            int matrixSize = diagMatrix.GetLength(0);
            int kCount = matrixSize * (matrixSize + 1) / 2;
            GH_Structure<GH_Integer> diagonals = new GH_Structure<GH_Integer>();
            for(int k = 0; k < kCount; k++)
            {
                //to access the UpperMatrix
                int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                int j = k + i * (i + 1) / 2 - matrixSize * i;
                if (diagMatrix[i, j] > 0)
                {
                    GH_Path current = new GH_Path(i);
                    GH_Path other = new GH_Path(j);
                    diagonals.Append(new GH_Integer(j), current);
                    diagonals.Append(new GH_Integer(i), other);
                }
            }
            return diagonals;
        }

            /// <summary>
            /// Creates a tree of the indices of the opposite corners of the valid diagonals of the polygon that connect each corner.
            /// </summary>
            /// <param name="diagMatrix">the matrix containing the diagonals of the polygon</param>
            /// <returns>Returns a tree of the indices of the opposite corner of each diagonal that connect to each corner.</returns>
            public DataTree<int> IndexTree(double[,] diagMatrix)
        {
            int matrixSize = diagMatrix.GetLength(0);
            int kCount = matrixSize * (matrixSize + 1) / 2;
            DataTree<int> diagonals = new DataTree<int>();
            for (int k = 0; k < kCount; k++)
            {
                //to access the UpperMatrix
                int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                int j = k + i * (i + 1) / 2 - matrixSize * i;
                if (diagMatrix[i, j] > 0)
                {
                    GH_Path current = new GH_Path(i);
                    GH_Path other = new GH_Path(j);
                    diagonals.Add(j, current);
                    diagonals.Add(i, other);
                }
            }
            return diagonals;
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.InternalDiagonals_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("a336b97c-c9fa-4c97-9d4a-750d88a4627b"); }
        }
    }
}
