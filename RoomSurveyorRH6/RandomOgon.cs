using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace RoomSurveyorRH6
{
    public class RandomOgon : GH_Component
    {
        public RandomOgon()
          : base("Random Ogon Generator", "Ogon1",
            "Random Ogon Generator - CutPaste. This algorithm generates a random orthogonal polygon by cutting or pasting rectangles to a seed rectangle. Is follows a similiar approach to one described by Tomas and Bajuelos 2004 InflatePaste and InflateCut algorithms, with some improvements to adapt it to floating point and reduce the likelyhood of failure in cuts. It implements one rule to prevent generating features that are smaller than a given dimension. This algorithm may fail if a point is outside the bounding box of the ogon or in certain interior very particular points.",
            "RoomSurveyor", "Polygon Generators")
        { 
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddRectangleParameter("Region", "R", "A rectangular region within which the ogon will be generated", GH_ParamAccess.item);
            pManager.AddNumberParameter("Minimum Feature","F","The length of the smallest feature", GH_ParamAccess.item, 0.10);
            pManager.AddBooleanParameter("Area Bias", "B", "If true the algorithm will favor Cuts or Pastes that minimize the area increase or decrease of the ogon", GH_ParamAccess.item, true);
            pManager.AddIntegerParameter("Number of sides", "N", "The number of sides of the generated Ogon. Only pair numbers are allowed", GH_ParamAccess.item, 6);
            pManager.AddIntegerParameter("Seed", "S", "An integer to be used as a seed for the random generator. If nothing is provided the system clock is used. Provide a seed if you wish to debug your algorithm", GH_ParamAccess.item);
            pManager[4].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Ogon", "O", "An orthogonal polygon with n sides", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.RandomOgon_Icon;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("8dd480ad-0f62-4770-9dd1-232f645fbd82"); }
        }
    }
}
