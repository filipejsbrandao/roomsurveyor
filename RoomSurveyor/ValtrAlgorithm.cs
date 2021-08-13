/*
 * MIT License
 * 
 * Copyright (c) 2017 Sander Verdonschot
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public static class ValtrAlgorithm
    {

        private static Random RAND = new Random();

        public static List<Point2d> GenerateRandomConvexPolygon(int n)
        {
            // Generate two lists of random X and Y coordinates
            List<Double> xPool = new List<Double>(n);
            List<Double> yPool = new List<Double>(n);

            for (int i = 0; i < n; i++)
            {
                xPool.Add(RAND.NextDouble());
                yPool.Add(RAND.NextDouble());
            }

            // Sort them
            xPool.Sort();
            yPool.Sort();

            // Isolate the extreme points
            Double minX = xPool[0];
            Double maxX = xPool[n - 1];
            Double minY = yPool[0];
            Double maxY = yPool[n - 1];

            // Divide the interior points into two chains & Extract the vector components
            List<Double> xVec = new List<Double>(n);
            List<Double> yVec = new List<Double>(n);

            double lastTop = minX, lastBot = minX;

            for (int i = 1; i < n - 1; i++)
            {
                double x1 = xPool[i];

                if (RAND.Next(2) == 0)
                {
                    xVec.Add(x1 - lastTop);
                    lastTop = x1;
                }
                else
                {
                    xVec.Add(lastBot - x1);
                    lastBot = x1;
                }
            }

            xVec.Add(maxX - lastTop);
            xVec.Add(lastBot - maxX);

            double lastLeft = minY, lastRight = minY;

            for (int i = 1; i < n - 1; i++)
            {
                double y1 = yPool[i];

                if (RAND.Next(2) == 0)
                {
                    yVec.Add(y1 - lastLeft);
                    lastLeft = y1;
                }
                else
                {
                    yVec.Add(lastRight - y1);
                    lastRight = y1;
                }
            }

            yVec.Add(maxY - lastLeft);
            yVec.Add(lastRight - maxY);

            // Randomly pair up the X- and Y-components
            yVec.Shuffle();
            //Collections.shuffle(yVec);

            // Combine the paired up components into vectors
            List<Point2d> vec = new List<Point2d>(n);

            for (int i = 0; i < n; i++)
            {
                vec.Add(new Point2d(xVec[i], yVec[i]));
            }

            // Sort the vectors by angle           
            vec.Sort((a, b) => (Math.Atan2(a.X, a.Y).CompareTo(Math.Atan2(b.X, b.Y))));

            //Collections.sort(vec, Comparator.comparingDouble(v->Math.atan2(v.getY(), v.getX())));

            // Lay them end-to-end
            double x = 0, y = 0;
            double minPolygonX = 0;
            double minPolygonY = 0;
            List<Point2d> points = new List<Point2d>(n);

            for (int i = 0; i < n; i++)
            {
                points.Add(new Point2d(x, y));

                x += vec[i].X;
                y += vec[i].Y;

                minPolygonX = Math.Min(minPolygonX, x);
                minPolygonY = Math.Min(minPolygonY, y);
            }

            // Move the polygon to the original min and max coordinates
            double xShift = minX - minPolygonX;
            double yShift = minY - minPolygonY;

            for (int i = 0; i < n; i++)
            {
                Point2d p = points[i];
                points[i] = new Point2d(p.X + xShift, p.Y + yShift);
            }

            return points;
        }
    }

    public static class ThreadSafeRandom
    {
        [ThreadStatic] private static Random Local;

        public static Random ThisThreadsRandom
        {
            get { return Local ?? (Local = new Random(unchecked(Environment.TickCount * 31 + Thread.CurrentThread.ManagedThreadId))); }
        }
    }

    static class MyExtensions
    {
        public static void Shuffle<T>(this IList<T> list)
        {
            int n = list.Count;
            while (n > 1)
            {
                n--;
                int k = ThreadSafeRandom.ThisThreadsRandom.Next(n + 1);
                T value = list[k];
                list[k] = list[n];
                list[n] = value;
            }
        }
    }
}
