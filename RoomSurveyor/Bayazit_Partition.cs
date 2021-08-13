using System;
using System.Collections;
using System.Collections.Generic;
using Rhino.Geometry;


/// <summary>
/// This class is adapted from Leonardo Temperaza implementation of Mark Bayazit decomposition algorithm as posted on the Unity forum
/// https://answers.unity.com/questions/977416/2d-polygon-convex-decomposition-code.html
/// I ported it to Rhino, replacing all Unity specifics with equivalent Rhino funtionality. All floats were replaced by doubles, the Vector2 class was replaced by Vector2d
/// for greater compatibility, but this comes with some limitations in the construction of objects when calling the class since Rhino Vector2d Constructor is very limited.
/// It should be noted that although the class outputs a List of Vector2d, these should be treated as 2D points and not 2D vectors.
/// Note to future-self: try to replace the Vector2d class with Point2d
/// Other more recent version is here: https://github.com/VelcroPhysics/VelcroPhysics/tree/master/VelcroPhysics/Tools
/// </summary>

/// <summary>
/// (By Leonardo Temperanza)
/// This class is taken from the "FarseerUnity" physics engine, which uses Mark Bayazit's decomposition algorithm.
/// https://github.com/pbhogan/FarseerUnity/blob/master/FarseerUnity/Farseer/Common/Decomposition/BayazitDecomposer.cs
/// I also have to make it work with self-intersecting polygons, so I'll use another different algorithm to decompose a self-intersecting polygon into several simple polygons,
/// and then I would decompose each of them into convex polygons.
/// </summary>

//From phed rev 36

/// <summary>
/// Convex decomposition algorithm created by Mark Bayazit (http://mnbayazit.com/)
/// For more information about this algorithm, see http://mnbayazit.com/406/bayazit
/// </summary>
namespace RoomSurveyor
{
    public static class Bayazit_Partition
    {

        private static Vector2d At(int i, List<Vector2d> vertices)
        {
            int s = vertices.Count;
            return vertices[i < 0 ? s - (-i % s) : i % s];
        }

        private static List<Vector2d> Copy(int i, int j, List<Vector2d> vertices)
        {
            List<Vector2d> p = new List<Vector2d>();
            while (j < i) j += vertices.Count;
            //p.reserve(j - i + 1);
            for (; i <= j; ++i)
            {
                p.Add(At(i, vertices));
            }
            return p;
        }

        /// <summary>
        /// Decompose the polygon into several smaller non-concave polygon.
        /// If the polygon is already convex, it will return the original polygon, unless it is over Settings.MaxPolygonVertices.
        /// Precondition: Counter Clockwise polygon
        /// </summary>
        /// <param name="vertices"></param>
        /// <returns></returns>
        public static List<List<Vector2d>> ConvexPartition(List<Vector2d> vertices)
        {
            //We force it to CCW as it is a precondition in this algorithm.
            ForceCounterClockWise(vertices);

            List<List<Vector2d>> list = new List<List<Vector2d>>();
            double d, lowerDist, upperDist;
            Vector2d p;
            Vector2d lowerInt = new Vector2d();
            Vector2d upperInt = new Vector2d(); // intersection points
            int lowerIndex = 0, upperIndex = 0;
            List<Vector2d> lowerPoly, upperPoly;

            for (int i = 0; i < vertices.Count; ++i)
            {
                if (Reflex(i, vertices))
                {
                    lowerDist = upperDist = double.MaxValue; // std::numeric_limits<qreal>::max();
                    for (int j = 0; j < vertices.Count; ++j)
                    {
                        // if line intersects with an edge
                        if (Left(At(i - 1, vertices), At(i, vertices), At(j, vertices)) &&
                            RightOn(At(i - 1, vertices), At(i, vertices), At(j - 1, vertices)))
                        {
                            // find the point of intersection
                            p = LineIntersect(At(i - 1, vertices), At(i, vertices), At(j, vertices),
                                                          At(j - 1, vertices));
                            if (Right(At(i + 1, vertices), At(i, vertices), p))
                            {
                                // make sure it's inside the poly
                                d = SquareDist(At(i, vertices), p);
                                if (d < lowerDist)
                                {
                                    // keep only the closest intersection
                                    lowerDist = d;
                                    lowerInt = p;
                                    lowerIndex = j;
                                }
                            }
                        }

                        if (Left(At(i + 1, vertices), At(i, vertices), At(j + 1, vertices)) &&
                            RightOn(At(i + 1, vertices), At(i, vertices), At(j, vertices)))
                        {
                            p = LineIntersect(At(i + 1, vertices), At(i, vertices), At(j, vertices),
                                                          At(j + 1, vertices));
                            if (Left(At(i - 1, vertices), At(i, vertices), p))
                            {
                                d = SquareDist(At(i, vertices), p);
                                if (d < upperDist)
                                {
                                    upperDist = d;
                                    upperIndex = j;
                                    upperInt = p;
                                }
                            }
                        }
                    }

                    // if there are no vertices to connect to, choose a point in the middle
                    if (lowerIndex == (upperIndex + 1) % vertices.Count)
                    {
                        Vector2d sp = ((lowerInt + upperInt) / 2);

                        lowerPoly = Copy(i, upperIndex, vertices);
                        lowerPoly.Add(sp);
                        upperPoly = Copy(lowerIndex, i, vertices);
                        upperPoly.Add(sp);
                    }
                    else
                    {
                        double highestScore = 0, bestIndex = lowerIndex;
                        while (upperIndex < lowerIndex) upperIndex += vertices.Count;
                        for (int j = lowerIndex; j <= upperIndex; ++j)
                        {
                            if (CanSee(i, j, vertices))
                            {
                                double score = 1 / (SquareDist(At(i, vertices), At(j, vertices)) + 1);
                                if (Reflex(j, vertices))
                                {
                                    if (RightOn(At(j - 1, vertices), At(j, vertices), At(i, vertices)) &&
                                        LeftOn(At(j + 1, vertices), At(j, vertices), At(i, vertices)))
                                    {
                                        score += 3;
                                    }
                                    else
                                    {
                                        score += 2;
                                    }
                                }
                                else
                                {
                                    score += 1;
                                }
                                if (score > highestScore)
                                {
                                    bestIndex = j;
                                    highestScore = score;
                                }
                            }
                        }
                        lowerPoly = Copy(i, (int)bestIndex, vertices);
                        upperPoly = Copy((int)bestIndex, i, vertices);
                    }
                    list.AddRange(ConvexPartition(lowerPoly));
                    list.AddRange(ConvexPartition(upperPoly));
                    return list;
                }
            }

            // polygon is already convex
            list.Add(vertices);

            //The polygons are not guaranteed to be without collinear points. We remove
            //them to be sure.
            for (int i = 0; i < list.Count; i++)
            {
                //list[i] = SimplifyTools.CollinearSimplify(list[i], 0);
            }

            //Remove empty vertice collections
            for (int i = list.Count - 1; i >= 0; i--)
            {
                if (list[i].Count == 0)
                    list.RemoveAt(i);
            }

            return list;
        }

        private static bool CanSee(int i, int j, List<Vector2d> vertices)
        {
            if (Reflex(i, vertices))
            {
                if (LeftOn(At(i, vertices), At(i - 1, vertices), At(j, vertices)) &&
                    RightOn(At(i, vertices), At(i + 1, vertices), At(j, vertices))) return false;
            }
            else
            {
                if (RightOn(At(i, vertices), At(i + 1, vertices), At(j, vertices)) ||
                    LeftOn(At(i, vertices), At(i - 1, vertices), At(j, vertices))) return false;
            }
            if (Reflex(j, vertices))
            {
                if (LeftOn(At(j, vertices), At(j - 1, vertices), At(i, vertices)) &&
                    RightOn(At(j, vertices), At(j + 1, vertices), At(i, vertices))) return false;
            }
            else
            {
                if (RightOn(At(j, vertices), At(j + 1, vertices), At(i, vertices)) ||
                    LeftOn(At(j, vertices), At(j - 1, vertices), At(i, vertices))) return false;
            }
            for (int k = 0; k < vertices.Count; ++k)
            {
                if ((k + 1) % vertices.Count == i || k == i || (k + 1) % vertices.Count == j || k == j)
                {
                    continue; // ignore incident edges
                }
                Vector2d intersectionPoint;
                if (LineIntersect2(At(i, vertices), At(j, vertices), At(k, vertices), At(k + 1, vertices), out intersectionPoint))
                {
                    return false;
                }
            }
            return true;
        }

        // precondition: ccw
        private static bool Reflex(int i, List<Vector2d> vertices)
        {
            return Right(i, vertices);
        }

        private static bool Right(int i, List<Vector2d> vertices)
        {
            return Right(At(i - 1, vertices), At(i, vertices), At(i + 1, vertices));
        }

        private static bool Left(Vector2d a, Vector2d b, Vector2d c)
        {
            return Area(ref a, ref b, ref c) > 0;
        }

        private static bool LeftOn(Vector2d a, Vector2d b, Vector2d c)
        {
            return Area(ref a, ref b, ref c) >= 0;
        }

        private static bool Right(Vector2d a, Vector2d b, Vector2d c)
        {
            return Area(ref a, ref b, ref c) < 0;
        }

        private static bool RightOn(Vector2d a, Vector2d b, Vector2d c)
        {
            return Area(ref a, ref b, ref c) <= 0;
        }

        private static double SquareDist(Vector2d a, Vector2d b)
        {
            double dx = b.X - a.X;
            double dy = b.Y - a.Y;
            return dx * dx + dy * dy;
        }

        //forces counter clock wise order.
        private static void ForceCounterClockWise(List<Vector2d> vertices)
        {
            if (!IsCounterClockWise(vertices))
            {
                vertices.Reverse();
            }
        }

        private static bool IsCounterClockWise(List<Vector2d> vertices)
        {
            //We just return true for lines
            if (vertices.Count < 3)
                return true;

            return (GetSignedArea(vertices) > 0.0);
        }

        //gets the signed area.
        private static double GetSignedArea(List<Vector2d> vertices)
        {
            int i;
            double area = 0;

            for (i = 0; i < vertices.Count; i++)
            {
                int j = (i + 1) % vertices.Count;
                area += vertices[i].X * vertices[j].Y;
                area -= vertices[i].Y * vertices[j].X;
            }
            area /= 2.0;
            return area;
        }

        //From Mark Bayazit's convex decomposition algorithm
        private static Vector2d LineIntersect(Vector2d p1, Vector2d p2, Vector2d q1, Vector2d q2)
        {
            Vector2d i = Vector2d.Zero;
            double a1 = p2.Y - p1.Y;
            double b1 = p1.X - p2.X;
            double c1 = a1 * p1.X + b1 * p1.Y;
            double a2 = q2.Y - q1.Y;
            double b2 = q1.X - q2.X;
            double c2 = a2 * q1.X + b2 * q1.Y;
            double det = a1 * b2 - a2 * b1;

            if (!DoubleEquals(det, 0))
            {
                // lines are not parallel
                i.X = (b2 * c1 - b1 * c2) / det;
                i.Y = (a1 * c2 - a2 * c1) / det;
            }
            return i;
        }

        //from Eric Jordan's convex decomposition library, it checks if the lines a0->a1 and b0->b1 cross.
        //if they do, intersectionPoint will be filled with the point of crossing. Grazing lines should not return true.
        private static bool LineIntersect2(Vector2d a0, Vector2d a1, Vector2d b0, Vector2d b1, out Vector2d intersectionPoint)
        {
            intersectionPoint = Vector2d.Zero;

            if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1)
                return false;

            double x1 = a0.X;
            double y1 = a0.Y;
            double x2 = a1.X;
            double y2 = a1.Y;
            double x3 = b0.X;
            double y3 = b0.Y;
            double x4 = b1.X;
            double y4 = b1.Y;

            //AABB early exit
            if (Math.Max(x1, x2) < Math.Min(x3, x4) || Math.Max(x3, x4) < Math.Min(x1, x2))
                return false;

            if (Math.Max(y1, y2) < Math.Min(y3, y4) || Math.Max(y3, y4) < Math.Min(y1, y2))
                return false;

            double ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3));
            double ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3));
            double denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
            if (Math.Abs(denom) < Rhino.RhinoMath.ZeroTolerance)
            {
                //Lines are too close to parallel to call
                return false;
            }
            ua /= denom;
            ub /= denom;

            if ((0 < ua) && (ua < 1) && (0 < ub) && (ub < 1))
            {
                intersectionPoint.X = (x1 + ua * (x2 - x1));
                intersectionPoint.Y = (y1 + ua * (y2 - y1));
                return true;
            }

            return false;
        }

        private static bool DoubleEquals(double value1, double value2)
        {
            return Math.Abs(value1 - value2) <= Rhino.RhinoMath.ZeroTolerance;
            //return Math.Abs(value1 - value2) <= Double.Epsilon;
        }

        // Returns a positive number if c is to the left of the line going from a to b. Positive number if point is left,
        //negative if point is right,and 0 if points are collinear.</returns>
        private static double Area(Vector2d a, Vector2d b, Vector2d c)
        {
            return Area(ref a, ref b, ref c);
        }

        //returns a positive number if c is to the left of the line going from a to b. Positive number if point is left, negative if point is right, and 0 if points are collinear.</returns>
        private static double Area(ref Vector2d a, ref Vector2d b, ref Vector2d c)
        {
            return a.X * (b.Y - c.Y) + b.X * (c.Y - a.Y) + c.X * (a.Y - b.Y);
        }

        //removes all collinear points on the polygon.
        private static List<Vector2d> CollinearSimplify(List<Vector2d> vertices, double collinearityTolerance)
        {
            //We can't simplify polygons under 3 vertices
            if (vertices.Count < 3)
                return vertices;

            List<Vector2d> simplified = new List<Vector2d>();

            for (int i = 0; i < vertices.Count; i++)
            {
                int prevId = PreviousIndex(vertices, i);
                int nextId = NextIndex(vertices, i);

                Vector2d prev = vertices[prevId];
                Vector2d current = vertices[i];
                Vector2d next = vertices[nextId];

                //If they collinear, continue
                if (Collinear(ref prev, ref current, ref next, collinearityTolerance))
                    continue;

                simplified.Add(current);
            }

            return simplified;
        }

        //gets the previous index.
        private static int PreviousIndex(List<Vector2d> vertices, int index)
        {
            if (index == 0)
            {
                return vertices.Count - 1;
            }
            return index - 1;
        }

        //nexts the index.
        private static int NextIndex(List<Vector2d> vertices, int index)
        {
            if (index == vertices.Count - 1)
            {
                return 0;
            }
            return index + 1;
        }

        private static bool Collinear(ref Vector2d a, ref Vector2d b, ref Vector2d c, double tolerance)
        {
            return DoubleInRange(Area(ref a, ref b, ref c), -tolerance, tolerance);
        }

        //checks if a floating point Value is within a specified range of values (inclusive).
        private static bool DoubleInRange(double value, double min, double max)
        {
            return (value >= min && value <= max);
        }
    }
    
}
