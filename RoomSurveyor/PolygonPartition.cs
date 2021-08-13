using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyor
{
    public class PolygonPartition : GH_Component
    {
        public PolygonPartition()
          : base("PolygonPartition", "PP",
            "Divide a non-convex polygon into convex components",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("P", "Polygon", "a simple non-intersecting polygon to partition", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("CP", "Convex Polygons", "A list of non-overlapping convex subdivisions", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Polyline poly = new Polyline();
            Curve curve = poly.ToNurbsCurve();
            List<Polyline> polygons = new List<Polyline>();

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

            List<Vector2d> vertices = new List<Vector2d>();

            for (int i = 0; i < poly.Count - 1; i++)
            {
                vertices.Add(new Vector2d(poly[i].X, poly[i].Y));
            }

            List<List<Vector2d>> convexPolygonPoints = Bayazit_Partition.ConvexPartition(vertices);

            foreach (List<Vector2d> polygonPts in convexPolygonPoints)
            {
                Polyline polygon = new Polyline();

                for (int i = 0; i < polygonPts.Count; i++)
                {
                    polygon.Add(new Point3d(polygonPts[i].X, polygonPts[i].Y, 0));
                    if (i == polygonPts.Count - 1)
                    {
                        polygon.Add(new Point3d(polygonPts[0].X, polygonPts[0].Y, 0));
                    }
                }
                polygons.Add(polygon);
            }

            DA.SetDataList(0, polygons);
        }

        /// <summary>
        /// Rebuilds a polyline with a list of <paramref name="vectors"/>
        /// </summary>
        /// <returns>The poly.</returns>
        /// <param name="polyline">Polyline.</param>
        /// <param name="vectors">Vectors.</param>
        private Polyline RebuildPoly(Polyline polyline, List<Vector2d> vectors)
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
                    Vector3d v = new Vector3d(vectors[j].X, vectors[j].Y, 0);
                    p = newPoly[j] + v;
                    newPoly.Add(p);
                }
            }
            return newPoly;
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("47802b12-490f-4e19-92aa-6f1036da22fd"); }
        }
    }
}
