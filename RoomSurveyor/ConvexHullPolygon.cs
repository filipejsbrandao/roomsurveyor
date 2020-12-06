using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino;

namespace RoomSurveyor
{
    public class ConvexHullPolygon : GH_Component
    {
        public override Grasshopper.Kernel.GH_Exposure Exposure
        {
            get { return GH_Exposure.hidden; }
        }

        public ConvexHullPolygon()
          : base("Convex Hull Simple Polygon", "CH SP",
            "Returns the convex hull (a convex polygon) of a simple polygon, i.e. any non-intersecting polygon. This component implements the O(n) algorithm described by Lee 1983. It is useful to obtain convex polygons of non-convex ones.",
            "RoomSurveyor", "Utils")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polygon", "P", "A closed simple planar polygon", GH_ParamAccess.item);
            pManager.AddPlaneParameter("Plane", "Pl", "The plane on which the polygon lies", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Convex Hull", "C", "A convex polygon", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Plane plane = Plane.Unset;
            Polyline poly = new Polyline();
            Polyline convexHull = new Polyline();
            Curve curve = poly.ToNurbsCurve();

            if (!DA.GetData(0, ref curve)) return;
            if (!DA.GetData(1, ref plane)) return;

            if (!curve.IsPlanar())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon must be a planar polyline");
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

            if (curvePlane.Normal.IsParallelTo(plane.Normal, 0.002) == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The polygon is not parallel to the provided plane");
                return;
            }

            //This will remove colinear lines bellow an angle threshold and will reposition the start point of the curve to a corner if that point is colinear.
            //It will also remove 0 length segments.
            poly.MergeColinearSegments(RhinoDoc.ActiveDoc.ModelAngleToleranceRadians / 100, true);

            poly.Transform(Transform.PlanarProjection(plane)); //We ensure it is on the plane
            poly = OrientPolygon.OrientPoly(poly, plane); //we make sure it is CCW oriented

            Transform transform = Transform.Identity;//We change nothing
            Transform reverseTrans = Transform.Identity;//We change nothing

            if (!plane.Equals(Plane.WorldXY))
            {
                transform = Transform.PlaneToPlane(plane, Plane.WorldXY);//otherwise we transform the polygon (place is on XY)
                reverseTrans = Transform.PlaneToPlane(Plane.WorldXY, plane);//And put it back in place when we are finished
            }

            poly.Transform(transform);

            convexHull.AddRange(ConvexHull.Lee83(poly));

            convexHull.Transform(reverseTrans);

            DA.SetData(0, convexHull);
        }

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
            get { return new Guid("cf378d32-98ef-4e72-a59b-7b4cdd6a68b3"); }
        }
    }
}
