using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyorRH6
{
    public class RoomSurveyStrict : GH_Component
    {

        public RoomSurveyStrict()
          : base("RoomSurveyStrict", "RS_S",
            "RoomSurvey version 7 - An algorithm to assist the survey of room plans using diagonals measured by the user between room corners. This version only requests the shortest diagonals",
            "RoomSurveyor", "RoomSurvey")
        {
        }
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A polyline that is garanteed to be a closed, non-intersecting and oriented polygon. Use the CheckPolygon component", GH_ParamAccess.item);
            pManager.AddNumberParameter("Side Lengths (m)", "SL", "A list containing the length in meters of each side in an anticlockwise sequence", GH_ParamAccess.list);
            pManager.AddNumberParameter("Diagonals (m)", "D", "A list containing the length of the diagonals in meters in the sequence requested by the Out of this component", GH_ParamAccess.list);
            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Output", "Out", "A list of string containing the sequence of needed diagonals to close the polygon", GH_ParamAccess.list);
            pManager.AddCurveParameter("Room", "R", "The surveyed room", GH_ParamAccess.item);
            pManager.AddLineParameter("Diagonals", "D", "Required Diagonals as lines for representation purposes", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            List<double> lengths = new List<double>();
            List<double> diagonals = new List<double>();

            if (!DA.GetData(0, ref curve)) return;
            if ((lengths != null) && !DA.GetDataList(1, lengths)) return;
            DA.GetDataList(2, diagonals);

            foreach (double l in lengths)
            {
                if (Math.Abs(l) < Double.Epsilon)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "All sides of the room must be larger than 0");
                    return;
                }
            }

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
            if (!poly.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "A closed polyline must be supplied");
                return;
            }
            if (poly.SegmentCount != lengths.Count)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "List of side lengths must have the same number of entries as the number of polyline sides");
                return;
            }

            List<Vector3d> polyVec = new List<Vector3d>();
            List<string> outText = new List<string>();
            List<Line> diagLines = new List<Line>();
            List<int> isTriVec = new List<int>();

            Polyline rebuiltPoly = new Polyline();
            double tol = 0.015; //in case the lenghts are scaled to mm this would became an int
            double error;

            poly = OrientPoly(poly);//Confirm that the polyline is CCW oriented

            //Contruct the vector chain that represents the new polyline
            for (int i = 0; i < poly.Count - 1; i++)
            {
                Vector3d v = new Vector3d(poly[i + 1] - poly[i]);
                v.Unitize();
                v *= lengths[i];
                polyVec.Add(v);
            }

            // verificar se a cadeia de vetores criará um poligono fechado
            if (IsClosed(polyVec, tol))
            {
                // GOT A SOLUTION
                error = ClosingError(polyVec) * 1000;
                rebuiltPoly = RebuildPoly(poly, polyVec);
                outText.Add("The Polygon is closed with a " + error + " mm error");
            }
            else
            {
                //LABEL THE POINTS ----------------------------------------------------------------------------------------------------------------------------------------------
                Rhino.RhinoDoc doc = Rhino.RhinoDoc.ActiveDoc;
                for (int i = 0; i < poly.Count - 1; i++)
                {
                    AddAnnotationText(doc, i.ToString(), poly[i], 1);
                }

                //SET THE VECTORS TO NOT FIXED -----------------------------------------------------------------------------------------------------------------------
                //We will use this to lock the corner's angles
                for (int j = 0; j < polyVec.Count; j++)
                    isTriVec.Add(0);

                //SET THE PROCESSING ORDER --------------------------------------------------------------------------------------------------------------------------------------

                //When a user provides a topology that already displays non-orthogonal angles is safe to assume that he is aware that some angles are non-orthogonal
                //It is also safe to assume that the provided non-ortho angle is not acurate, since the user does not have tools to measure non-ortho angles on site.
                //Surveying interior spaces with hand measuring tools requires in certain instances two people to complete the task. With laser distance meters this is not the case
                //but if the user only has flexible or rigid tape measures tools longer measurements might require two users. 
                //We know that longest diagonals reduce the impact of measurement errors for laser distance meters.
                //On the other hand, for flexible tape measures it becomes more difficult to properly stretch the tape with larger distances.
                //It is always preferable to lay the flexible tape on the ground, but this may not be possible for diagonal measurements when there are obstacles.
                //Metalic tape measures, are less susceptible to the previous problem, but frequently have smaller lenghts (2, 3, 5, less frequently 8m long)
                //Consequently, if the user needs to measure a distance larger that the tape's length, he will need to make several measurements and add up the distances.
                //The easiest way to do this is to lay the tape on the floor, mark the end of the measurement and continue. The problem with this approach is that 
                //it can introduce several errors: 
                //1. As the measurement is being partitioned the user cannot use the tape as a visual guide to garantee that the measurement is being taken along a straight line.
                //2. The added calculation step can also introduce user errors.
                //As we can't make assumptions on which tool is being used to measure the room we should favor measuring non ortho angles first. 
                //As we are aiming for lead users (DIYers, small contractors or architects) we shall assume that they have adequate tools and know how to use them. 
                //Furthermore, laser distance meters are widely available, inexpensive tools and have become a very common tool in practice. 

                //-----1 the non-orthogonal angles
                //-----2 the longest diagonals

                //lets create a vector with the order of the points
                List<int> orderedPoints = new List<int>();
                orderedPoints.AddRange(OrderByAngle(poly));

                //Then the function that creates a matrix with all the diagonals of the polygon
                double[,] diagMatrix = UpperMatrix(poly);

                List<int> orderedDiagonals = new List<int>();
                //if (IsPolyOrtho(poly))
                //{
                //    orderedDiagonals.AddRange(DiagonalOrder(diagMatrix));
                //}
                //else
                //{
                    orderedDiagonals.AddRange(DiagonalOrder(diagMatrix, orderedPoints));
                    orderedDiagonals.AddRange(DiagonalOrder(diagMatrix));
                    orderedDiagonals = orderedDiagonals.Distinct().ToList();
                //}

                //SUGEST THE FIRST DIAGONALS --------------------------------------------------------------------------------------------------------------------------------------
                for (int c = 0; c < orderedDiagonals.Count; c++)
                {//create a diagonal using the diagonal vector for the user to see and a text intruction indicating the action
                    int k = orderedDiagonals[c];
                    int matrixSize = diagMatrix.GetLength(0);
                    int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                    int j = k + i * (i + 1) / 2 - matrixSize * i;
                    Line diagonal = new Line(poly[i], poly[j]);
                    diagLines.Add(diagonal);
                    outText.Add("Measure the distance from Point " + i + " to Point " + j);
                    if (diagonals.Count > c)
                    {//if the user has supplied a diagonal add it to the LowerMatrix
                        diagMatrix[j, i] = diagonals[c];
                        bool ijIsTri = false;
                        bool jiIsTri = false;
                        //Vector3d ijChain = new Vector3d();
                        //Vector3d jiChain = new Vector3d();

                        /*
                        for (int l = i; l < j; l++) { ijChain += polyVec[l]; }
                        for (int n = j; n < polyVec.Count; n++) { jiChain += polyVec[n]; }
                        for (int m = 0; m < i; m++) { jiChain += polyVec[m]; }
                        // A tolerance of 7,5mm causes deviations of around 2.2 cm in a 6.175m diagonal
                        // most recent Laser Distance Meters are rated with a +/- 1.5mm accuracy, in our case user introduced errors are likely higher
                        */
                        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //--------INTERNAL TRIANGLES---------------------------------------------------
                        //the triangles formed by two edges of the polygon and a diagonal
                        //the triangles formed by one edge of the polygon and two diagonals
                        //the internal triangles of the polygon formed by three diagonals
                        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if (InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out bool itoj) >= 0 && diagonals[c] > 0)
                        {
                            int P1 = InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out itoj);
                            int P0 = (itoj) ? i : j;
                            int P2 = (itoj) ? j : i;
                            double isleft = IsLeft(poly[P0], poly[P1], poly[P2]);
                            if (!Triangulate(P0, P1, P2, itoj, diagonals[c], isTriVec, polyVec, isleft))
                            {
                                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The length of the last provided diagonal " + diagonals[c] + " is longer than the sum of the length of its adjacent walls from " + P0 + " to " + P2);
                                return;
                            }
                            if (itoj) { ijIsTri = true; } else { jiIsTri = true; }
                        }

                        //check if the already provided diagonals will triangulate something
                        int storedDiagonal = StoredDiagonal(diagMatrix, isTriVec, out itoj);
                        if (storedDiagonal > 0)
                        {
                            while (storedDiagonal > 0 && Count0(isTriVec.ToArray(), matrixSize - 1, matrixSize, matrixSize) > 1)
                            {
                                i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * storedDiagonal + 0.25));
                                j = storedDiagonal + i * (i + 1) / 2 - matrixSize * i;
                                int Pone = InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out itoj);
                                int Pzero = (itoj) ? i : j;
                                int Ptwo = (itoj) ? j : i;
                                double isL = IsLeft(poly[Pzero], poly[Pone], poly[Ptwo]);
                                if (!Triangulate(Pzero, Pone, Ptwo, itoj, diagMatrix[j, i], isTriVec, polyVec, isL))
                                {
                                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The length of the last provided diagonal " + diagMatrix[j, i] + " is longer than the sum of the length of its adjacent walls from " + Pzero + " to " + Ptwo);
                                    return;
                                }
                                if (itoj) { ijIsTri = true; } else { jiIsTri = true; }
                                storedDiagonal = StoredDiagonal(diagMatrix, isTriVec, out itoj);
                            }
                        }
                        //test if the diagonal[c] closes the poligonal chain i to j or the polygonal chain j to i
                        //Now we need to check if any or both of the polygonal chains became closed with the diagonal
                        /*if (Math.Abs(ijChain.Length - diagonals[c]) < tol / 2 && diagonals[c] > 0)
                        {
                            ijIsTri = true;//lock the vectors from i+1 to j
                            for (int n = i + 1; n < j; n++)
                            {
                                isTriVec[n] = 1;
                            }
                        }

                        if (Math.Abs(jiChain.Length - diagonals[c]) < tol / 2 && diagonals[c] > 0)
                        {
                            jiIsTri = true;//lock the vectors of the polygonal chain from j+1 to i
                            for (int n = j + 1; n < polyVec.Count; n++)
                            {
                                isTriVec[n] = 1;
                            }

                            for (int m = 0; m < i; m++)
                            {
                                isTriVec[m] = 1;
                            }
                        }*/

                        if (ijIsTri && !jiIsTri || !ijIsTri && jiIsTri || diagonals[c] < 0 || !ijIsTri && !jiIsTri)
                        {
                            //remove the diagonals between triangulated points from the diagonalOrder list
                            if (ijIsTri && !jiIsTri || !ijIsTri && jiIsTri || diagonals[c] < 0)
                            {
                                RemoveDiagonals(ijIsTri, jiIsTri, i, j, diagMatrix, orderedDiagonals);
                            }
                            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            //--------101 PATTERN-------------
                            //101 paterns can speedup the triangulation and help to triangulate concave corners. The aim is always to fix the 0 in the sequence.
                            //we first create an ordered list of valid diagonals from large to small. 
                            //Request the user the largest valid diagonal and if not possible then request one of the remaining ones.
                            //if no valid diagonals exist then continue.
                            //There can be multiple 101 patterns, the key is to find the first zero on both sides: 0...101.0
                            //if the center 0 is i we need to find the i-nth 0 and the i+nth 0 and all the possible diagonals inside.
                            //order the diagonals by size
                            //CURRENTLY DIAGONALORDER101 ONLY RETURNS 4 DIAGONALS
                            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            int pattern101 = Sequence101(isTriVec.ToArray(), isTriVec.Count);//HERE ADD CHECKS FOR -1 value
                            if (pattern101 >= 0)
                            {
                                int center = pattern101;
                                List<int> diagonals101 = new List<int>();
                                diagonals101.AddRange(DiagonalOrder101(center, diagMatrix));
                                k = StoredDiagonal(diagMatrix, diagonals101);//returns the k index of the stored diagonal on the lower matrix or -1 if none found
                                if (k < 0)
                                {//determine if the first diagonal is on the list (but not before the current c) and remove it
                                    if (diagonals101.Any())
                                    {
                                        int found = orderedDiagonals.IndexOf(diagonals101[0], c + 1);
                                        if (found != -1)
                                        {
                                            orderedDiagonals.RemoveAt(found);
                                        }//insert it on the next c
                                        orderedDiagonals.Insert(c + 1, diagonals101[0]);
                                    }
                                }
                            }
                            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            //OOO Pattern - search for a 000 unique sequence; find the two intersections of the last not fixed edges and select the correct one
                            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            int P0 = Pattern000(isTriVec.ToArray(), isTriVec.Count, out int P2, out int P1);
                            if (P0 >= 0)
                            {
                                Vector3d p0p1 = new Vector3d(0, 0, 0);
                                Vector3d p1p2 = new Vector3d(0, 0, 0);
                                Vector3d p2p0 = new Vector3d(0, 0, 0);
                                for (int l = P0; l < P1; l++) { p0p1 += polyVec[l]; }
                                for (int n = P1; n < P2; n++) { p1p2 += polyVec[n]; }
                                for (int m = P2; m < polyVec.Count; m++) { p2p0 += polyVec[m]; }
                                for (int p = 0; p < P0; p++) { p2p0 += polyVec[p]; }
                                if (!IsValidTriangle(p0p1.Length, p1p2.Length, p2p0.Length, out double diff))
                                {
                                    goto Skip;
                                }
                                Vector3d newP1P2 = NextVector(Point3d.Origin, p0p1, p0p1.Length, p1p2.Length, p2p0.Length);
                                if (IsLeft(poly[P0], poly[P1], poly[P2]) < 0)
                                {
                                    Vector3d per = p0p1;
                                    per.Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                                    newP1P2.Transform(Transform.Mirror(Point3d.Origin, per));
                                }
                                newP1P2.Reverse();
                                double angle1 = AngleBetweenVectors(p1p2, newP1P2, false);
                                for (int n = P1; n < P2; n++)
                                {
                                    Vector3d v = polyVec[n];
                                    v.Transform(Transform.Rotation(angle1, Point3d.Origin));
                                    polyVec[n] = v;
                                }
                                newP1P2.Reverse();
                                p1p2 = newP1P2;
                                Vector3d newP2P0 = NextVector(Point3d.Origin, p1p2, p1p2.Length, p2p0.Length, p0p1.Length);
                                if (IsLeft(poly[P1], poly[P2], poly[P0]) < 0)
                                {
                                    Vector3d per = p1p2;
                                    per.Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                                    newP2P0.Transform(Transform.Mirror(Point3d.Origin, per));
                                }
                                newP2P0.Reverse();
                                angle1 = AngleBetweenVectors(p2p0, newP2P0, false);
                                for (int n = P2; n < polyVec.Count; n++)
                                {
                                    Vector3d v = polyVec[n];
                                    v.Transform(Transform.Rotation(angle1, Point3d.Origin));
                                    polyVec[n] = v;
                                }
                                for (int m = 0; m < P0; m++)
                                {
                                    Vector3d v = polyVec[m];
                                    v.Transform(Transform.Rotation(angle1, Point3d.Origin));
                                    polyVec[m] = v;
                                }
                                isTriVec[P0] = 1;
                                isTriVec[P1] = 1;
                                isTriVec[P2] = 1;
                                error = ClosingError(polyVec) * 1000;
                                rebuiltPoly = RebuildPoly(poly, polyVec);
                                outText.Add("The Polygon is closed with a " + error + " mm error");
                                break;
                            }
                        }
                    Skip:
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //OO PATTERN - search for a 00 unique sequence; find the vector from the start to the end of the pattern
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        int startP00 = Pattern00(isTriVec.ToArray(), isTriVec.Count, out int endP00);
                        if (startP00 >= 0)
                        {
                            Vector3d ij = new Vector3d(0, 0, 0);
                            Vector3d ji = new Vector3d(0, 0, 0);
                            for (int n = startP00; n < endP00; n++) { ij += polyVec[n]; }
                            for (int m = endP00; m < polyVec.Count; m++) { ji += polyVec[m]; }
                            for (int l = 0; l < startP00; l++) { ji += polyVec[l]; }
                            double angle = AngleBetweenVectors(ij, ji, false);
                            for (int n = startP00; n < endP00; n++)
                            {
                                Vector3d v = polyVec[n];
                                v.Transform(Transform.Rotation(angle, Point3d.Origin));
                                polyVec[n] = v;
                            }
                            error = ClosingError(polyVec) * 1000;
                            rebuiltPoly = RebuildPoly(poly, polyVec);
                            outText.Add("The Polygon is closed with a " + error + " mm error");
                            break;
                        }
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        //If all corners are triangulated, end processing
                        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if (Count0(isTriVec.ToArray(), matrixSize - 1, matrixSize, matrixSize) <= 1)
                        {
                            error = ClosingError(polyVec) * 1000;
                            rebuiltPoly = RebuildPoly(poly, polyVec);
                            outText.Add("The Polygon is closed with a " + error + " mm error");
                            break;
                        }
                    }
                    else
                    {
                        //Debug Method
                        rebuiltPoly.Add(poly[0]);
                        Point3d pt = new Point3d();

                        for (int m = 0; m < polyVec.Count; m++)
                        {
                            pt = rebuiltPoly[m] + polyVec[m];
                            rebuiltPoly.Add(pt);
                        }
                        //End Degug method
                        break;
                    }
                }
            }
            DA.SetDataList(0, outText);
            DA.SetData(1, rebuiltPoly);
            DA.SetDataList(2, diagLines);

        }
        /// <summary>
        /// Triangulate the specified polygonal chain polyVec between P0 and P2, with center on P1 and diagonal as the length of the diagonal.
        /// </summary>
        /// <returns>The triangulate.</returns>
        /// <param name="P0">The start point of the diagonal.</param>
        /// <param name="P1">The corner to triangulate.</param>
        /// <param name="P2">The end point of the diagonal.</param>
        /// <param name="itoj">If set to <c>true</c> the triangulation occurs on the ij polygonal chain of the polygon.</param>
        /// <param name="diagonal">the length of the diagonal.</param>
        /// <param name="isTriVec">A list of <list type="int"></list> containing the state of triangulation of each corner of the polygon.</param>
        /// <param name="polyVec">A list of the chain of vectors of the polygon.</param>
        public static bool Triangulate(int P0, int P1, int P2, bool itoj, double diagonal, List<int> isTriVec, List<Vector3d> polyVec, double isleft)
        {
            int matrixSize = polyVec.Count;
            Vector3d p0p1 = new Vector3d(0, 0, 0);
            Vector3d p1p2 = new Vector3d(0, 0, 0);
            if (itoj)
            {
                for (int n = P0; n < P1; n++) { p0p1 += polyVec[n]; }
                for (int m = P1; m < P2; m++) { p1p2 += polyVec[m]; }
            }
            else
            {
                if (P1 < P0)
                {
                    for (int n = P0; n < matrixSize; n++) { p0p1 += polyVec[n]; }
                    for (int m = 0; m < P1; m++) { p0p1 += polyVec[m]; }
                    for (int l = P1; l < P2; l++) { p1p2 += polyVec[l]; }
                }
                else
                {
                    for (int n = P0; n < P1; n++) { p0p1 += polyVec[n]; }
                    for (int m = P1; m < matrixSize; m++) { p1p2 += polyVec[m]; }
                    for (int l = 0; l < P2; l++) { p1p2 += polyVec[l]; }
                }
            }
            Vector3d newp1p2 = NextVector(Point3d.Origin, p0p1, p0p1.Length, p1p2.Length, diagonal);
            if (!newp1p2.IsValid)
            {
                return false;
            }
            if (isleft < 0)
            {
                Vector3d per = p0p1;
                per.Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                newp1p2.Transform(Transform.Mirror(Point3d.Origin, per));
            }
            newp1p2.Reverse();
            double newAngle = AngleBetweenVectors(p1p2, newp1p2, false);
            int lastFixed = Sequence1N_CCW(isTriVec.ToArray(), P1, isTriVec.Count, 0);
            if (lastFixed >= 0)
            {
                RotateFixedChain(lastFixed, P1, polyVec, newAngle);
            }
            else
            {
                RotateFixedChain(P1, P1, polyVec, newAngle);
            }
            isTriVec[P1] = 1;

            return true;
        }

        /// <summary>
        /// Ises the valid triangle.
        /// </summary>
        /// <returns><c>true</c>, if valid triangle was ised, <c>false</c> otherwise.</returns>
        /// <param name="length1">Length1.</param>
        /// <param name="length2">Length2.</param>
        /// <param name="length3">Length3.</param>
        /// <param name="error">Error.</param>
        public static bool IsValidTriangle(double length1, double length2, double length3, out double error)
        {
            bool triangle = true;
            error = -1;

            double[] sides = { length1, length2, length3 };
            Array.Sort(sides);
            if (sides[0] + sides[1] < sides[2])
            {
                error = sides[2] - (sides[0] + sides[1]);
                return false;
            }

            return triangle;
        }

        /// <summary>
        /// Removes the diagonals.
        /// </summary>
        /// <param name="diagonalMatrix">Diagonal matrix.</param>
        /// <param name="orderedDiagonals">Ordered diagonals.</param>
        /// <param name="isTriVec">Is tri vec.</param>
        public static void ReorderDiagonals(double[,] diagonalMatrix, List<int> orderedDiagonals, List<int> isTriVec)
        {
            int k;
            int matrixSize = diagonalMatrix.GetLength(0);

            //What is a diagonal that does not contribute to triangulate any corner?
            //Answer: A diagonal that goes from i to j such that all corners between i and j are triangulated corners.
            //Proof: If there is a diagonal from i to j, such that i, j and all corners from i to j are triangulated, that can be used to triangulate the chain from j to i, then any of the diagonals that have triangulated the corners from i to j may also be used to triangulate j to i.
            //Are there any other cases?
            //Yes. If i and j are non-triangulated corners and all corners between i and j are triangulated, the diagonal will not contribute to triangulate any corner between i and j.


            //++++++++++++TODO+++++++++++++++++++++++++
            //A better method to remove diagonals by looking at patterns of 0s and 1s

        }
        /// <summary>
        /// Removes the diagonals between fixed corners of the ij or ji polygonal chain.
        /// </summary>
        /// <param name="ijIsTriangulated">If set to <c>true</c> ij is triangulated.</param>
        /// <param name="jiIsTriangulated">If set to <c>true</c> ji is triangulated.</param>
        /// <param name="i">The index of the from point of the diagonal.</param>
        /// <param name="j">The index of the to point of the diagonal.</param>
        /// <param name="diagonalMatrix">The diagonal matrix.</param>
        /// <param name="orderedDiagonals">The list of ordered diagonals.</param>
        public static void RemoveDiagonals(bool ijIsTriangulated, bool jiIsTriangulated, int i, int j, double[,] diagonalMatrix, List<int> orderedDiagonals)
        {
            int matrixSize = diagonalMatrix.GetLength(0);

            int k;
            if (ijIsTriangulated)
            {
                for (int n = i; n < j; n++)
                {
                    for (int m = n + 2; m <= j; m++)
                    {
                        if (!(n == i && m == j) && diagonalMatrix[m, n] == 0)
                        {
                            k = n * matrixSize - n * (n + 1) / 2 + m;
                            orderedDiagonals.Remove(k);
                        }
                    }
                }
            }
            else if (jiIsTriangulated)
            {
                for (int n = j; n < matrixSize; n++)
                {
                    for (int m = n + 2; m < matrixSize; m++)
                    {
                        if (!(n == i && m == j) && diagonalMatrix[m, n] == 0) //here we try to avoid removing diagonals that the user already provided
                        {
                            k = n * matrixSize - n * (n + 1) / 2 + m;
                            orderedDiagonals.Remove(k);
                        }
                    }
                }
                for (int n = 0; n <= i; n++)
                {
                    for (int m = n + 2; m <= i; m++)
                    {
                        if (diagonalMatrix[m, n] <= 0) //here we try to avoid removing diagonals that the user already provided
                        {
                            k = n * matrixSize - n * (n + 1) / 2 + m;
                            orderedDiagonals.Remove(k);
                        }
                    }
                    for (int m = j; m < matrixSize; m++)
                    {
                        if (!(n == i && m == j) && diagonalMatrix[m, n] == 0) //here we try to avoid removing diagonals that the user already provided
                        {
                            k = n * matrixSize - n * (n + 1) / 2 + m;
                            orderedDiagonals.Remove(k);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Searches for user provided diagonals in the LowerMatrix
        /// </summary>
        /// <returns>Returns the first that is found</returns>
        /// <param name="diagonalMatrix">Diagonal matrix.</param>
        /// <param name="diagonals101">The list of diagonals to search (ordered) </param>
        public static int StoredDiagonal(double[,] diagonalMatrix, List<int> diagonals101)
        {
            int found = -1;
            int n = diagonals101.Count;
            int m = 0;
            int matrixSize = diagonalMatrix.GetLength(0);

            while (m < n)
            {
                int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * diagonals101[m] + 0.25));
                int j = diagonals101[m] + i * (i + 1) / 2 - matrixSize * i;
                if (diagonalMatrix[j, i] > 0)
                {
                    found = diagonals101[m];
                    return found;
                }
                m++;
            }

            return found;
        }

        /// <summary>
        /// Search for a stored diagonal that can triangulate the polygon and return the k index.
        /// </summary>
        /// <returns>The k index of the diagonal on the upper matrix.</returns>
        /// <param name="diagonalMatrix">Diagonal matrix.</param>
        /// <param name="isTriVec">The list of triangulation state of the polygon vectors.</param>
        /// <param name="ij">If set to <c>true</c> the diagonal triangulates the polygon on the ij chain.</param>
        public static int StoredDiagonal(double[,] diagonalMatrix, List<int> isTriVec, out bool ij)
        {
            int found = -1;
            ij = true;

            int matrixSize = diagonalMatrix.GetLength(0);
            int kCount = matrixSize * (matrixSize + 1) / 2;
            for (int k = 0; k < kCount; k++)
            {
                //to access the UpperMatrix
                int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                int j = k + i * (i + 1) / 2 - matrixSize * i;
                if (diagonalMatrix[j, i] > 0)
                {
                    if (InternalTriangle(isTriVec.ToArray(), i, j, matrixSize, out ij) >= 0)
                    {
                        return k;
                    }
                }
            }

            return found;
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //                         +++++++++++++++++++++++  PATTERNS +++++++++++++++++++++++
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// <summary>
        /// Test if a diagonal can triangulate the IJ chain or the JI chain of a polygon.
        /// </summary>
        /// <returns>The index of the corner that can be triangulated or -1 if none</returns>
        /// <param name="v">An array with 1s and 0s containing the state of triangulation of the polygon's corners</param>
        /// <param name="i">The start point.</param>
        /// <param name="j">The end point.</param>
        /// <param name="n">The length of the v array.</param>
        /// <param name="ij">If set to <c>true</c> the triangulation can be done on the ij chain.</param>
        public static int InternalTriangle(int[] v, int i, int j, int n, out bool ij)
        {
            ij = true;
            int P1 = -1;

            int count = 0;
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ IJ
            count = Count0(v, i, j, n);
            if (count == 1)
            {
                for (int m = i + 1; m < j; m++)
                {
                    if (v[m] == 0)
                    {
                        return m;
                    }
                }
            }
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ JI
            count = Count0(v, j, i, n);
            if (count == 1)
            {
                int start = (j + 1) % n;
                int end = (i == 0) ? n - 1 : i - 1;

                if (start <= end)
                {
                    for (int l = start; l <= end; l++)
                    {
                        if (v[l] == 0)
                        {
                            ij = false;
                            return l;
                        }
                    }
                }
                else
                {
                    for (int m = start; m < n; m++)
                    {
                        if (v[m] == 0)
                        {
                            ij = false;
                            return m;
                        }
                    }
                    for (int k = 0; k <= end; k++)
                    {
                        if (v[k] == 0)
                        {
                            ij = false;
                            return k;
                        }
                    }
                }

            }

            return P1;
        }

        /// <summary>
        /// Count the number of 0 in the specified v, P0, P2 and n.
        /// </summary>
        /// <returns>The number of 0s in the list.</returns>
        /// <param name="v">V.</param>
        /// <param name="P0">P0.</param>
        /// <param name="P2">P2.</param>
        /// <param name="n">N.</param>
        public static int Count0(int[] v, int P0, int P2, int n)
        {
            int count = 0;
            int start = (P0 + 1) % n;
            int end = (P2 == 0) ? n - 1 : P2 - 1;

            if (start <= end)
            {
                for (int i = start; i <= end; i++)
                {
                    if (v[i] == 0)
                        count++;
                }
            }
            else
            {
                for (int j = start; j < n; j++)
                {
                    if (v[j] == 0)
                        count++;
                }
                for (int k = 0; k <= end; k++)
                {
                    if (v[k] == 0)
                        count++;
                }
            }

            return count;
        }

        /// <summary>
        /// Finds the last sequence of 101 values in a list.
        /// </summary>
        /// <returns>The middle index of the pattern (the 0)</returns>
        /// <param name="v">An array of integers 1 or 0</param>
        /// <param name="n">The number of values in the array</param>
        public static int Sequence101(int[] v, int n)
        {
            int start = -1;
            if (n < 4)
            {
                return start;
            }
            for (int i = 0; i < n; i++)
            {
                int r = (i + 1) % n;
                int s = (i + 2) % n;
                if (v[i] == 1 && v[r] == 0 && v[s] == 1)
                {
                    start = r;
                }
            }
            return start;
        }

        /// <summary>
        /// Searches for the last index of a sequence of values that is different from a breaking_mark
        /// Useful for searching the end of a 101 pattern in the CCW direction or checking if the next corner is fixed. 
        /// It returns -1 if the next corner is not fixed (if it is equal to 0), it returns -1 if all
        /// the corners are not fixed (if all the values are equal to the breaking mark) or if all are fixed
        /// (if all corners are different from the breaking mark).
        /// </summary>
        /// <returns>The index of the last 1 or fixed angle or vector.</returns>
        /// <param name="v">An array of values 1 and 0. </param>
        /// <param name="P0">The center of the 101 pattern.</param>
        /// <param name="n">The number of elements in the array.</param>
        /// <param name="breaking_mark">0.</param>
        public static int Sequence1N_CCW(int[] v, int P0, int n, int breaking_mark)
        {

            int end = -1;
            int Pstart = (P0 + 1) % n;
            int Pend = (P0 == 0) ? n - 1 : P0 - 1;
            bool found = false;

            if (v[Pstart] == breaking_mark) { return end; }

            for (int i = Pstart; i < n; i++)
            {
                if (v[i] == breaking_mark)
                {
                    found = true;
                    return (i + n - 1) % n;
                }
            }
            if (!found)
            {
                for (int i = 0; i <= Pend; i++)
                {
                    if (v[i] == breaking_mark)
                    {
                        found = true;
                        return (i + n - 1) % n;
                    }
                }
            }
            return end;
        }

        /// <summary>
        /// Find any 000 Pattern.
        /// </summary>
        /// <returns>The index of the start, middle and end of the 000 pattern.</returns>
        /// <param name="v">V.</param>
        /// <param name="n">N.</param>
        /// <param name="end">End.</param>
        public static int Pattern000(int[] v, int n, out int end, out int middle)
        {
            int start = -1;
            end = -1;
            middle = -1;

            if (n < 3)
            {
                middle = -1;
                end = -1;
                return -1;
            }
            else if (n == 3)
            {
                if (v[0] == 0 && v[1] == 0 && v[2] == 0)
                {
                    middle = 1;
                    end = 2;
                    return 0;
                }
            }

            int zeros = 0;

            for (int i = 0; i < n; i++)
            {
                if (v[i] == 0)
                {
                    zeros++;
                    switch (zeros)
                    {
                        case 1:
                            start = i;
                            break;
                        case 2:
                            middle = i;
                            break;
                        case 3:
                            end = i;
                            break;
                        case 4:
                            middle = -1;
                            end = -1;
                            return -1;
                    }
                }
            }
            if (middle == -1 || end == -1)
            {
                middle = -1;
                end = -1;
                return -1;
            }

            return start;
        }

        /// <summary>
        /// Find any 00 Pattern.
        /// </summary>
        /// <returns>The index of the start and end of the 00 pattern.</returns>
        /// <param name="v">V.</param>
        /// <param name="n">N.</param>
        /// <param name="end">End.</param>
        public static int Pattern00(int[] v, int n, out int end)
        {
            int start = -1;
            end = -1;

            if (n < 3)
            {
                end = -1;
                return -1;
            }
            int zeros = 0;

            for (int i = 0; i < n; i++)
            {
                if (v[i] == 0)
                {
                    zeros++;
                    switch (zeros)
                    {
                        case 1:
                            start = i;
                            break;
                        case 2:
                            end = i;
                            break;
                        case 3:
                            end = -1;
                            return -1;
                    }
                }
            }
            if (end == -1)
            {
                return -1;
            }

            return start;
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        /// <summary>
        /// Tests if a diagonal is a valid diagonal of a polygon.
        /// </summary>
        /// <returns><c>true</c>, if valid diagonal was ised, <c>false</c> otherwise.</returns>
        /// <param name="diagonal">Diagonal.</param>
        /// <param name="polyline">Polyline.</param>
        public static bool IsValidDiagonal(Line diagonal, Polyline polyline)
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
            bool valid = false;

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

                double first = AngleAtCorner(current, previous, next); //must be Left so larger than 0
                double second = AngleAtCorner(current, diagonalEnd, next); //must be Left so larger than 0

                if (first > second)
                    valid = true;
            }

            return valid;
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        /// <summary>
        /// A method to find the set of valid diagonals that triangulate a 101 pattern in a closed polygon
        /// It takes the center of the pattern and the upper matrix containing all the diagonals of the polygon.
        /// </summary>
        /// <returns>Returns an ordered list of valid diagonals, from large to small, that triangulate a 101 pattern in a closed polygon.</returns>
        /// <param name="center">The center of the 101 pattern, or the index of the corner to be triangulated</param>
        /// <param name="diagMatrix">The upper diagonal matrix containing all the valid diagonals of the polygon</param>
        public static List<int> DiagonalOrder101(int center, double[,] diagMatrix)
        {

            List<int> diagonalOrder = new List<int>();
            List<KeyValuePair<int, double>> diagonalVector = new List<KeyValuePair<int, double>>();
            int n = diagMatrix.GetLength(0);

            int start = center - 2;
            if (center == 0)
            {
                start = n - 2;
            }
            else if (center == 1)
            {
                start = n - 1;
            }

            int next = (start + 1) % n;
            int next3 = (start + 3) % n;
            int next4 = (start + 4) % n;

            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    if ((i == start || i == next || i == next3 || i == next4) && (j == start || j == next || j == next3 || j == next4))
                    {
                        if (diagMatrix[i, j] > 0 && diagMatrix[j, i] >= 0) //here we are checking if the diagonal is a valid diagonal of the polygon and if the user has not rejected the diagonal
                        {
                            int k = i * n - i * (i + 1) / 2 + j;
                            KeyValuePair<int, double> diag = new KeyValuePair<int, double>(k, diagMatrix[i, j]);
                            diagonalVector.Add(diag);
                        }
                    }
                }
            }

            diagonalVector.Sort(delegate (KeyValuePair<int, double> firstPair, KeyValuePair<int, double> nextPair)
            {
                return nextPair.Value.CompareTo(firstPair.Value);
            }
              );

            diagonalOrder = (from d in diagonalVector select d.Key).ToList();


            return diagonalOrder;
        }

        /// <summary>
        /// Creates a list of integers k, ordered by diagonal length from large to small, of the valid diagonals of the polygon.
        /// The index k can be used to access the [i,j] of the matrix
        /// </summary>
        /// <param name="diagMatrix">the matrix containing the diagonals of the polygon</param>
        /// <returns>Returns a list of indices with the order of the diagonals (large to smallest).</returns>
        public static List<int> DiagonalOrder(double[,] diagMatrix)
        {
            List<int> diagonalOrder = new List<int>();
            int matrixSize = diagMatrix.GetLength(0);
            int kCount = matrixSize * (matrixSize + 1) / 2;
            List<KeyValuePair<int, double>> diagonalVector = new List<KeyValuePair<int, double>>();
            for (int k = 0; k < kCount; k++)
            {
                //to access the UpperMatrix
                int i = (int)Math.Floor(matrixSize + 0.5 - Math.Sqrt(matrixSize * (matrixSize + 1) - 2 * k + 0.25));
                int j = k + i * (i + 1) / 2 - matrixSize * i;
                KeyValuePair<int, double> diag = new KeyValuePair<int, double>(k, diagMatrix[i, j]);
                if (diag.Value > 0)
                {
                    diagonalVector.Add(diag);
                }
            }

            diagonalVector.Sort(delegate (KeyValuePair<int, double> firstPair, KeyValuePair<int, double> nextPair)
            {
                return nextPair.Value.CompareTo(firstPair.Value);
            }
              );

            diagonalOrder = (from d in diagonalVector select d.Key).ToList();

            return diagonalOrder;
        }

        /// <summary>
        /// Return the diagonals from points i-1 to i+1 for each point i that is not ortho and whose angle is smaller than Pi.
        /// </summary>
        /// <param name="diagMatrix"></param>
        /// <param name="angleOrder"></param>
        /// <returns></returns>
        public static List<int> DiagonalOrder(double[,] diagMatrix, List<int> angleOrder)
        {
            List<int> diagonalOrder = new List<int>();
            int matrixSize = diagMatrix.GetLength(0);

            for (int n = 0; n < angleOrder.Count; n++)
            {
                int current = angleOrder[n];
                int i = current - 1;
                int j = current + 1;
                if (current == matrixSize - 1)
                {
                    j = i;
                    i = 0;
                }
                else if (current == 0)
                {
                    i = j;
                    j = matrixSize - 1;
                }
                if (diagMatrix[i, j] > 0.0) //should I return -1 if there is no valid diagonal?
                {// calculate k for i,j
                    int k = i * matrixSize - i * (i + 1) / 2 + j;
                    diagonalOrder.Add(k);
                }
            }

            return diagonalOrder;
        }
        /// <summary>
        /// Fills the upper part of a diagonal matrix with the valid diagonals of a closed polyline
        /// </summary>
        /// <param name="poly"> A valid closed CCW oriented polyline </param>
        /// <returns>A matrix with the values of the diagonals of a polygon. 
        /// Diagonals to the current one, the next and the previous point return 0.
        /// Invalid diagonals return -1 </returns>
        public static double[,] UpperMatrix(Polyline poly)
        {
            double[,] diagMatrix = new double[poly.Count - 1, poly.Count - 1];

            int n = poly.Count - 1;

            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    if (i == 0)
                    {
                        if (j == i || j == i + 1 || j == poly.Count - 2)
                        {
                            diagMatrix[i, j] = 0;
                        }
                        else
                        {
                            Line diagonal = new Line(poly[i], poly[j]);
                            if (IsValidDiagonal(diagonal, poly))
                            {
                                diagMatrix[i, j] = diagonal.Length;
                            }
                            else
                            {
                                diagMatrix[i, j] = -1;
                            }
                        }
                    }
                    else
                    {
                        if (j == i || j == i + 1)
                        {
                            diagMatrix[i, j] = 0;
                        }
                        else
                        {
                            Line diagonal = new Line(poly[i], poly[j]);
                            if (IsValidDiagonal(diagonal, poly))
                            {
                                diagMatrix[i, j] = diagonal.Length;
                            }
                            else
                            {
                                diagMatrix[i, j] = -1;
                            }
                        }
                    }
                }
            }

            return diagMatrix;
        }
        /// <summary>
        /// Sets the processing order of the points of the polyline by angle
        /// The angles are sorted from small to large
        /// </summary>
        /// <returns>A list of indices in which to process the corners of the polyline.</returns>
        /// <param name="poly">A closed polyline that represents the room to triangulate.</param>
        public static List<int> OrderByAngle(Polyline poly)  //I have also created a method that does the angle check using the vectors instead of the points
        {
            List<int> pointOrder = new List<int>();
            List<KeyValuePair<int, double>> angles_pts = new List<KeyValuePair<int, double>>();

            for (int i = 0; i < poly.Count - 1; i++)
            {
                if (i == 0)
                {
                    angles_pts.Add(new KeyValuePair<int, double>(i, AngleAtCorner(poly[i], poly[poly.Count - 2], poly[i + 1])));
                }
                else
                {
                    angles_pts.Add(new KeyValuePair<int, double>(i, AngleAtCorner(poly[i], poly[i - 1], poly[i + 1])));
                }
            }

            angles_pts.Sort(delegate (KeyValuePair<int, double> firstPair, KeyValuePair<int, double> nextPair)
            {
                return firstPair.Value.CompareTo(nextPair.Value);
            }
                    );

            List<KeyValuePair<int, double>> orthoAngles = new List<KeyValuePair<int, double>>();
            for (int j = 0; j < angles_pts.Count; j++)
            {
                if ((Math.Abs(angles_pts[j].Value - Math.PI / 2) < 0.0001) || (Math.Abs(angles_pts[j].Value - 3 * Math.PI / 2) < 0.0001))
                {
                    orthoAngles.Add(angles_pts[j]);
                    angles_pts.RemoveAt(j);
                    j -= 1;
                }
            }
            angles_pts.AddRange(orthoAngles); //if it is needed to add ortho angles to the list...

            pointOrder = (from a in angles_pts select a.Key).ToList();

            return pointOrder;
        }

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// <summary>
        /// Calculates the CCW angle at a corner given the previous and the next point
        /// </summary>
        /// <returns>Returns the angle in radians</returns>
        /// <param name="currPt">current corner - the corner</param>
        /// <param name="prevPt">Previous corner - the previous point or the other end of the edge P1 to P0</param>
        /// <param name="nextPt">Next corner - the next point or the other end of the edge P0 tp P2</param>
        private static double AngleAtCorner(Point3d currPt, Point3d prevPt, Point3d nextPt)
        {
            double a = currPt.DistanceTo(prevPt);
            double b = currPt.DistanceTo(nextPt);
            double c = prevPt.DistanceTo(nextPt);

            double angle = FindAngleSSS(a, b, c);
            double isLeft = IsLeft(prevPt, nextPt, currPt);

            if (isLeft > 0)
            { // is a concave angle
                angle = 2 * Math.PI - angle;
            }
            else if (Math.Abs(isLeft) < double.Epsilon)
            { // is colinear
                angle = Math.PI;
            }
            return angle;
        }

        /// <summary>
        /// Assuming that the two vectors supplied are the translation vectors in a CCW oriented polygon
        /// i.e.: the first vector currVect is the translation vector from i to i+1 and the second vector 
        /// prevVect is the translation vector i-1 to i, we want to know the internal angle of the polygon at the corner i.
        /// </summary>
        /// <returns>The <see cref="T:System.Double"/>.</returns>
        /// <param name="currVect">The current vector.</param>
        /// <param name="prevVect">The previous vector.</param>
        /// <param name="returndegrees">If set to <c>true</c> return degrees else returns radians.</param>
        public static double AngleBetweenVectors(Vector3d currVect, Vector3d prevVect, bool returndegrees)
        {
            double toppart = 0;

            Point3d P0 = new Point3d(prevVect);
            Point3d P1 = new Point3d(currVect);
            Point3d P2 = Point3d.Origin; //a point with (0,0,0)coordinates

            prevVect.Reverse();//we need to reverse the previous vector because the vectors are CCW oriented

            for (int d = 0; d < 3; d++) toppart += currVect[d] * prevVect[d];

            double currVect2 = 0; //currVect squared
            double prevVect2 = 0; //prevVect squared
            for (int d = 0; d < 3; d++)
            {
                currVect2 += currVect[d] * currVect[d];
                prevVect2 += prevVect[d] * prevVect[d];
            }

            double bottompart = 0;
            bottompart = Math.Sqrt(currVect2 * prevVect2);
            double division = (toppart / bottompart > 1) ? 1 : toppart / bottompart;

            double angle = Math.Acos(division);

            if (IsLeft(P0, P1, P2) < 0)
            {
                angle = 2 * Math.PI - angle;
            }
            if (returndegrees) angle *= 360.0 / (2 * Math.PI);
            return angle;
        }

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// <summary>
        /// Given a closed polyline this method returns a Counter Clockwise oriented polyline
        /// </summary>
        /// <param name="poly"> A closed polyline</param>
        /// <returns>A CCW oriented polyline</returns>
        public static Polyline OrientPoly(Polyline poly)
        {

            Vector3d unitZ = new Vector3d(0, 0, 1);
            Curve a = poly.ToNurbsCurve();
            if (a.ClosedCurveOrientation(unitZ) == CurveOrientation.Clockwise)
                poly.Reverse();

            return poly;
        }

        /// <summary>
        /// Returns the tolerance required for closing the polygon
        /// </summary>
        /// <returns>The effective tolerance in mm.</returns>
        /// <param name="vectors">A List of vectors that create the polygon.</param>
        private double ClosingError(List<Vector3d> vectors)
        {
            Vector3d s = new Vector3d();
            foreach (Vector3d v in vectors) s += v;
            double error = s.Length;
            return error;
        }

        /// <summary>
        /// Rebuilds a polyline with a list of <paramref name="vectors"/>
        /// </summary>
        /// <returns>The poly.</returns>
        /// <param name="polyline">Polyline.</param>
        /// <param name="vectors">Vectors.</param>
        private Polyline RebuildPoly(Polyline polyline, List<Vector3d> vectors)
        {
            Polyline newPoly = new Polyline();
            newPoly.Add(polyline[0]);
            Point3d p = new Point3d();
            for (int j = 0; j < vectors.Count; j++)
            {
                if (j == vectors.Count - 1)
                {
                    p = newPoly[0];
                    newPoly.Add(p);
                }
                else
                {
                    p = newPoly[j] + vectors[j];
                    newPoly.Add(p);
                }
            }
            return newPoly;
        }

        /// <summary>
        /// Rotates the fixed CCW chain of vectors on the spot.
        /// </summary>
        /// <param name="lastFixed">The index of the last fixed vector or corner.</param>
        /// <param name="center">The index of the vector that was changed or the center of the 101 pattern.</param>
        /// <param name="polyVec">The list of vectors to be changed.</param>
        /// <param name="angle">The angle by which all the vectors need to be rotated.</param>
        public static void RotateFixedChain(int lastFixed, int center, List<Vector3d> polyVec, double angle)
        {

            if (lastFixed >= center)
            {
                for (int n = center; n <= lastFixed; n++)
                {
                    Vector3d v = polyVec[n];
                    v.Transform(Transform.Rotation(angle, Point3d.Origin));
                    polyVec[n] = v;
                }
            }
            else
            {
                for (int n = center; n < polyVec.Count; n++)
                {
                    Vector3d v = polyVec[n];
                    v.Transform(Transform.Rotation(angle, Point3d.Origin));
                    polyVec[n] = v;
                }
                for (int m = 0; m <= lastFixed; m++)
                {
                    Vector3d v = polyVec[m];
                    v.Transform(Transform.Rotation(angle, Point3d.Origin));
                    polyVec[m] = v;
                }
            }
        }


        /// <summary>
        /// Calculates the direction of the next polygon edge given the vector of the current edge, its origin, the lengths of both edges and the diagonal between
        /// the start point of the first edge and the end point of the second edge.
        /// </summary>
        /// <returns>The vector.</returns>
        /// <param name="startPt">Start point of the first edge. It can be the origin.</param>
        /// <param name="vector">The vector from the start point to the end point of the first edge.</param>
        /// <param name="firstLength">The lenght of the first side.</param>
        /// <param name="secondLength">The length of the second side.</param>
        /// <param name="diagonal">Diagonal between the origin and the end point of the second edge.</param>
        public static Vector3d NextVector(Point3d startPt, Vector3d vector, double firstLength, double secondLength, double diagonal)
        {
            Point3d nextPt = startPt + vector;
            double angle = FindAngleSSS(firstLength, secondLength, diagonal);
            Vector3d nextVec = vector;
            nextVec.Unitize();
            nextVec *= secondLength;
            nextVec.Transform(Transform.Rotation(Math.PI - angle, nextPt));

            return nextVec;
        }

        /// <summary>
        /// Finds the internal angle of the line with length <paramref name="a"/> with line with length <paramref name="b"/> given the length of the diagonal <paramref name="c"/>.
        /// </summary>
        /// <returns>The internal angle in radians</returns>
        /// <param name="a">The length of first side .</param>
        /// <param name="b">The length of the second side component.</param>
        /// <param name="c">The length of the diagonal</param>
        private static double FindAngleSSS(double a, double b, double c)
        {
            return Math.Acos((a * a + b * b - c * c) / (2 * a * b));
        }

        /// <summary>
        /// Determines if a third point is Left, On or Right of an infinite 2D line defined by two points
        /// </summary>
        /// <returns>Returns a positive number if it is Left, a negative number if it is Right and 0 if it is on the line</returns>
        /// <param name="P0">P0 - the first point</param>
        /// <param name="P1">P1 - the second point</param>
        /// <param name="P2">P2 - the point to test</param>
        private static double IsLeft(Point3d P0, Point3d P1, Point3d P2)
        {
            return ((P1.X - P0.X) * (P2.Y - P0.Y) - (P2.X - P0.X) * (P1.Y - P0.Y));
        }

        /// <summary>
        /// Checks if a list of vectors creates a closed polygon within a given tolerance.
        /// </summary>
        /// <returns><c>true</c>, if the resultant vector is smaller than the given <paramref name="tolerance"/>, <c>false</c> otherwise.</returns>
        /// <param name="vectors">A list contaning the vectors.</param>
        /// <param name="tolerance">Tolerance.</param>
        private bool IsClosed(List<Vector3d> vectors, double tolerance)
        {
            Vector3d s = new Vector3d();
            foreach (Vector3d v in vectors) s += v;
            //IsTiny method takes a double
            if (s.IsTiny(tolerance) == true || s.IsZero)
                return true;
            else
                return false;

        }
        /// <summary>
        /// Checks if all the angles in polyline are orthogonal, either 90 or 270.
        /// </summary>
        /// <returns><c>true</c>, if an ortho-polygon was input, <c>false</c> otherwise.</returns>
        /// <param name="polyline">Polyline.</param>
        private bool IsPolyOrtho(Polyline polyline)
        {
            double turn = 0.0;
            bool ortho = false;

            for (int i = 0; i < polyline.SegmentCount; i++)
            {
                Point3d p0, p1, p2;
                if (i == 0)
                {
                    p0 = polyline[polyline.Count - 2];
                    p1 = polyline[i];
                    p2 = polyline[i + 1];
                }
                else
                {
                    p0 = polyline[i - 1];
                    p1 = polyline[i];
                    p2 = polyline[i + 1];
                }
                Vector3d v1 = p1 - p0;
                Vector3d v2 = p2 - p1;
                turn += Math.Abs(v1 * v2);
            }

            if (Math.Abs(turn) < 0.00001) { ortho = true; }

            return ortho;
        }


        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //This method will need to be improved
        public static void AddAnnotationText(Rhino.RhinoDoc doc, string text, Point3d pt, double height)
        {
            Point3d iPt = pt;
            string iText = text;
            double iHeight = height;
            const string font = "Arial";
            Plane plane = doc.Views.ActiveView.ActiveViewport.ConstructionPlane();
            plane.Origin = iPt;
            Guid id = doc.Objects.AddText(iText, plane, iHeight, font, false, false);
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.RoomSurveyStrict_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("e6800037-f3a9-408b-8bc6-1e51d2236c0c"); }
        }
    }
}
