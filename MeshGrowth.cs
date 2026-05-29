using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using Plankton;
using PlanktonGh;
using KangarooSolver;
using KangarooSolver.Goals;
using System.Threading.Tasks;
using System.Linq;

/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
 private void RunScript(
    int xnum, int ynum,
    List<double> data,
    UVInterval interval2d,
    Mesh inputMesh,
    Mesh collisionMesh,
    double splitRatio,
    double remeshDecayRate,
    List<double> targetLengths,
    List<double> growthRates,
    double growthDecay,
    List<bool> isAnchored,
    double hingeStrength,
    bool reset,
    int viewNum,
    ref object A, ref object B, ref object C, ref object D, ref object E, ref object F)
{
    // Pre-processing 
    bool remeshComplete;
    double currentHingeStrength = hingeStrength;
    double springStrength       = 1.0;
    double selfPushDist1        = 1.0;
    double selfPushDist2        = 1.5;

    // Clamp splitRatio to [0, 1]
    double clampedSplitRatio = splitRatio >= 0 ? (splitRatio <= 1 ? splitRatio : 1) : 0;
    
    //  Reset & Initialise 
    if (reset)
    {
        // Build 2D environment data array from flat input list
        envData = new double[xnum][];
        for (int i = 0; i < xnum; i++)
        {
            envData[i] = new double[ynum];
            for (int j = 0; j < ynum; j++)
            {
                envData[i][j] = data[j + i * ynum];
            }
        }

        xCount  = xnum;
        yCount  = ynum;
        uvDomain = interval2d;
        collideMesh = collisionMesh;

        perVertexTargetLength = targetLengths;
        halfedgeMesh = inputMesh.ToPlanktonMesh();
        perVertexGrowthRate = growthRates;

        physicsSystem = new KangarooSolver.PhysicalSystem();
        anchorGoals   = new List<IGoal>();
        springGoals   = new List<IGoal>();
        hingeGoals    = new List<IGoal>();
        isVertexFixed = new List<bool>();
        _iterationCount = 0;

        int vertexIndex = 0;
        foreach (var vertex in halfedgeMesh.Vertices)
        {
            Point3d position    = vertex.ToPoint3d();
            double envValue     = Interpolate(position);
            perVertexEnvValue.Add(envValue);

            physicsSystem.AddParticle(position, 10);

            bool shouldBeFixed = isAnchored[vertexIndex] || envValue <= 0;
            if (!shouldBeFixed)
            {
                isVertexFixed.Add(false);
            }
            else
            {
                IGoal anchor = new AnchorXYZ(position, true, true, true, 100);
                anchor.PIndex = new int[1] { vertexIndex };
                anchorGoals.Add(anchor);
                isVertexFixed.Add(true);
            }

            vertexIndex++;
        }

        // Initialise spring and hinge goals for each edge
        int edgeCount = halfedgeMesh.Halfedges.Count / 2;
        for (int i = 0; i < edgeCount; i++)
        {
            var  startHalfedge = halfedgeMesh.Halfedges[2 * i];
            var  pairHalfedge  = halfedgeMesh.Halfedges[2 * i + 1];
            int  startVert     = startHalfedge.StartVertex;
            int  endVert       = pairHalfedge.StartVertex;

            double restLength = (halfedgeMesh.Vertices[startVert].ToPoint3d()
                               - halfedgeMesh.Vertices[endVert].ToPoint3d()).Length;
            springGoals.Add(new Spring(startVert, endVert, restLength, springStrength));

            double restAngle;
            int[] hingeVerts = GetHingeVertices(2 * i, halfedgeMesh, out restAngle);

            bool hasBothFaces = startHalfedge.AdjacentFace != -1 && pairHalfedge.AdjacentFace != -1;
            double hingeStr   = hasBothFaces ? currentHingeStrength : 0;
            hingeGoals.Add(new Hinge(hingeVerts[0], hingeVerts[1], hingeVerts[2], hingeVerts[3], 0, hingeStr));
        }
    }


    // Main Loop (Remesh + Grow) 
    double remeshTolerance = 1.0;
    int    remeshIterations = 0;

    double clampedRemeshDecay = remeshDecayRate;
    clampedRemeshDecay = clampedRemeshDecay < 0.96 ? (clampedRemeshDecay < 0 ? 0 : clampedRemeshDecay) : 0.96;

    double clampedGrowthDecay = growthDecay;
    clampedGrowthDecay = clampedGrowthDecay < 1 ? (clampedGrowthDecay < 0.6 ? 0.6 : clampedGrowthDecay) : 1;


    if (!reset)
    {
        //Remeshing Loop
        do
        {
            remeshComplete = true;

            var halfedges  = halfedgeMesh.Halfedges;
            int numEdges   = halfedges.Count / 2;
            var edgeLengths = halfedges.GetLengths();
            var vertices   = halfedgeMesh.Vertices;

            // Split too long
            for (int i = 0; i < numEdges; i++)
            {
                if (halfedges[2 * i].IsUnused) continue;

                int startVert = halfedges[2 * i].StartVertex;
                int endVert   = halfedges[2 * i + 1].StartVertex;

                // Only process edges where at least one endpoint is free
                if (isVertexFixed[startVert] && isVertexFixed[endVert]) continue;

                double avgTargetLength = (perVertexTargetLength[startVert] + perVertexTargetLength[endVert]) / 2.0;

                if (edgeLengths[2 * i] > (1 + remeshTolerance) * 1.25f * avgTargetLength)
                {
                    int newHalfedge = halfedges.TriangleSplitEdge(2 * i);
                    if (newHalfedge == -1) continue; // split failed

                    double growthStart = perVertexGrowthRate[startVert] *= clampedGrowthDecay;
                    double growthEnd   = perVertexGrowthRate[endVert]   *= clampedGrowthDecay;
                    double minGrowth   = growthStart > growthEnd ? growthEnd : growthStart;
                    double avgGrowth   = (growthStart + growthEnd) / 2.0;
                    double newGrowth   = (1 - clampedSplitRatio) * minGrowth + clampedSplitRatio * avgGrowth;

                    remeshComplete = false;

                    isVertexFixed.Add(false);
                    perVertexTargetLength.Add(avgTargetLength);
                    perVertexGrowthRate.Add(newGrowth);

                    // Place the new midpoint particle
                    Point3d startPt = vertices[startVert].ToPoint3d();
                    Point3d endPt   = vertices[endVert].ToPoint3d();
                    Point3d midPt   = (startPt + endPt) / 2.0;
                    physicsSystem.AddParticle(midPt, 10);
                    perVertexEnvValue.Add(Interpolate(midPt));

                    // Resize goal lists to match new halfedge count
                    int sizeDiff = halfedges.Count / 2 - springGoals.Count;
                    springGoals.AddRange(new Spring[sizeDiff]);
                    hingeGoals.AddRange(new Hinge[sizeDiff]);

                    // Rebuild goals for all edges around the new vertex
                    IEnumerable<int> neighborHalfedges = halfedges.GetVertexCirculator(newHalfedge);
                    D = halfedges[newHalfedge].StartVertex;

                    foreach (var neighborHe in neighborHalfedges)
                    {
                        int edgeIndex = neighborHe / 2;
                        double restAngle;
                        int[] hingeVerts = GetHingeVertices(2 * edgeIndex, halfedgeMesh, out restAngle);

                        bool hasBothFaces = halfedgeMesh.Halfedges[neighborHe].AdjacentFace != -1
                                         && halfedgeMesh.Halfedges[halfedgeMesh.Halfedges.GetPairHalfedge(neighborHe)].AdjacentFace != -1;
                        double hingeStr = hasBothFaces ? currentHingeStrength : 0;

                        hingeGoals[edgeIndex]  = new Hinge(hingeVerts[0], hingeVerts[1], hingeVerts[2], hingeVerts[3], 0, hingeStr);

                        double edgeLength = vertices[hingeVerts[0]].ToPoint3d().DistanceTo(vertices[hingeVerts[1]].ToPoint3d());
                        springGoals[edgeIndex] = new Spring(hingeVerts[0], hingeVerts[1], edgeLength, springStrength);
                    }
                }
            }

            halfedgeMesh.Compact();

            // Refresh references after compaction
            halfedges   = halfedgeMesh.Halfedges;
            numEdges    = halfedges.Count / 2;
            edgeLengths = halfedges.GetLengths();
            vertices    = halfedgeMesh.Vertices;

            // Collapse too short
            for (int i = 0; i < numEdges; i++)
            {
                if (halfedges[2 * i].IsUnused) continue;

                int startVert = halfedges[2 * i].StartVertex;
                int endVert   = halfedges[2 * i + 1].StartVertex;

                if (isVertexFixed[startVert] && isVertexFixed[endVert]) continue;

                double avgTargetLength = (perVertexTargetLength[startVert] + perVertexTargetLength[endVert]) / 2.0;

                if (edgeLengths[2 * i] < (1 - remeshTolerance) * 0.25f * avgTargetLength)
                {
                    int nextHe    = halfedges[2 * i].NextHalfedge;
                    int prevHe    = halfedges[2 * i + 1].PrevHalfedge;
                    int leftFace  = halfedges[2 * i].AdjacentFace;
                    int rightFace = halfedges[2 * i + 1].AdjacentFace;

                    int collapseStart = halfedges[2 * i].StartVertex;
                    int collapseEnd   = halfedges[2 * i + 1].StartVertex;

                    float midX = (vertices[collapseStart].X + vertices[collapseEnd].X) / 2f;
                    float midY = (vertices[collapseStart].Y + vertices[collapseEnd].Y) / 2f;
                    float midZ = (vertices[collapseStart].Z + vertices[collapseEnd].Z) / 2f;

                    double growthStart = perVertexGrowthRate[collapseStart];
                    double growthEnd   = perVertexGrowthRate[collapseEnd];
                    double keepGrowth  = growthStart > growthEnd ? growthStart : growthEnd;

                    int survivingHalfedge = halfedges.CollapseEdge(2 * i);
                    if (survivingHalfedge == -1) continue; // collapse failed

                    remeshComplete = false;

                    vertices.SetVertex(collapseStart, midX, midY, midZ);
                    perVertexTargetLength[collapseStart] = avgTargetLength;
                    perVertexGrowthRate[collapseStart]   = keepGrowth;
                    perVertexTargetLength.RemoveAt(collapseEnd);
                    perVertexGrowthRate.RemoveAt(collapseEnd);

                    // Remove goals for the collapsed edge and its faces' edges
                    springGoals.RemoveAt(i);
                    hingeGoals.RemoveAt(i);
                    physicsSystem.DeleteParticle(collapseEnd);

                    if (leftFace  != -1) { springGoals.RemoveAt(nextHe / 2); hingeGoals.RemoveAt(nextHe / 2); }
                    if (rightFace != -1) { springGoals.RemoveAt(prevHe / 2); hingeGoals.RemoveAt(prevHe / 2); }

                    // Remap vertex indices from removed vertex to surviving vertex
                    IEnumerable<int> neighborHalfedges = halfedges.GetVertexCirculator(survivingHalfedge);
                    foreach (var neighborHe in neighborHalfedges)
                    {
                        int edgeIndex = neighborHe / 2;

                        for (int k = 0; k < 4; k++)
                            if (hingeGoals[edgeIndex].PIndex[k] == collapseEnd)
                                hingeGoals[edgeIndex].PIndex[k] = collapseStart;

                        for (int k = 0; k < 2; k++)
                            if (springGoals[edgeIndex].PIndex[k] == collapseEnd)
                                springGoals[edgeIndex].PIndex[k] = collapseStart;

                        // Flag the opposite vertex as fixed to avoid processing stale edges
                        isVertexFixed[halfedges[halfedges.GetPairHalfedge(neighborHe)].StartVertex] = true;
                    }
                }
            }

            halfedgeMesh.Compact();

            // Flip
            int finalEdgeCount = halfedges.Count / 2;
            for (int i = 0; i < finalEdgeCount; i++)
            {
                if (halfedges[2 * i].IsUnused) continue;

                // Edge must be interior and all surrounding edges must also be interior
                bool isInteriorEdge =
                    halfedges[2 * i].AdjacentFace     != -1 &&
                    halfedges[2 * i + 1].AdjacentFace != -1 &&
                    halfedges[halfedges.GetPairHalfedge(halfedges[2 * i].NextHalfedge)].AdjacentFace     != -1 &&
                    halfedges[halfedges.GetPairHalfedge(halfedges[2 * i].PrevHalfedge)].AdjacentFace     != -1 &&
                    halfedges[halfedges.GetPairHalfedge(halfedges[2 * i + 1].NextHalfedge)].AdjacentFace != -1 &&
                    halfedges[halfedges.GetPairHalfedge(halfedges[2 * i + 1].PrevHalfedge)].AdjacentFace != -1;

                if (!isInteriorEdge) continue;

                int v1 = halfedges[2 * i].StartVertex;
                int v2 = halfedges[2 * i + 1].StartVertex;
                int v3 = halfedges[halfedges[halfedges[2 * i].NextHalfedge].NextHalfedge].StartVertex;
                int v4 = halfedges[halfedges[halfedges[2 * i + 1].NextHalfedge].NextHalfedge].StartVertex;

                Point3d P1 = vertices[v1].ToPoint3d();
                Point3d P2 = vertices[v2].ToPoint3d();
                Point3d P3 = vertices[v3].ToPoint3d();
                Point3d P4 = vertices[v4].ToPoint3d();

                // Sum of opposite angles for current diagonal
                double currentDiagonalAngles =
                    Vector3d.VectorAngle(new Vector3d(P3 - P1), new Vector3d(P4 - P1)) +
                    Vector3d.VectorAngle(new Vector3d(P4 - P2), new Vector3d(P3 - P2));

                // Sum of opposite angles for flipped diagonal
                double flippedDiagonalAngles =
                    Vector3d.VectorAngle(new Vector3d(P1 - P4), new Vector3d(P2 - P4)) +
                    Vector3d.VectorAngle(new Vector3d(P2 - P3), new Vector3d(P1 - P3));

                if (flippedDiagonalAngles > currentDiagonalAngles)
                {
                    remeshComplete = false;
                    halfedgeMesh.Halfedges.FlipEdge(2 * i);

                    springGoals[i] = new Spring(v3, v4, P3.DistanceTo(P4), springStrength);

                    double unusedAngle;
                    GetHingeVertices(2 * i, halfedgeMesh, out unusedAngle);
                    hingeGoals[i] = new Hinge(v3, v4, v2, v1, 0, currentHingeStrength);
                }
            }

            // Continue until converged or tolerance has decayed to near zero
            remeshComplete = remeshTolerance > 0.05 ? false : remeshComplete;
            remeshTolerance *= clampedRemeshDecay;
            remeshIterations++;

        } while (!remeshComplete && remeshIterations < 101);

        Print("{0}", remeshIterations - 1);
        halfedgeMesh.Compact();


        // Growth Step 

        var allGoals = new List<IGoal>();

        for (int i = 0; i < halfedgeMesh.Halfedges.Count / 2; i++)
        {
            Spring spring = springGoals[i] as Spring;
            int    sv     = spring.PIndex[0];
            int    ev     = spring.PIndex[1];

            // Grow rest length by combined growth and environment rates
            spring.RestLength += (perVertexGrowthRate[sv] + perVertexGrowthRate[ev])
                               * (perVertexEnvValue[sv]   + perVertexEnvValue[ev]) / 4.0;
            allGoals.Add(spring);

            Hinge hinge = hingeGoals[i] as Hinge;
            if (hinge.Strength != 0) allGoals.Add(hinge);
        }

        allGoals.AddRange(anchorGoals);
        allGoals.Add(new GrowthPush(selfPushDist1, selfPushDist2, halfedgeMesh));
        allGoals.Add(new MeshPush(halfedgeMesh, collideMesh, 0.5));

        physicsSystem.Step(allGoals, true, 0.1);

        Point3d[] updatedPositions = physicsSystem.GetPositionArray();

        for (int n = 0; n < updatedPositions.Count(); n++)
        {
            Vector3d displacement = updatedPositions[n] - halfedgeMesh.Vertices[n].ToPoint3d();
            perVertexEnvValue[n]  = Interpolate(updatedPositions[n]);

            // Freeze vertex if growth has effectively stopped or it left the growth field
            bool growthStopped  = !isVertexFixed[n] && displacement.Length < 0.0001 && perVertexGrowthRate[n] < 0.0001;
            bool outsideGrowthField = perVertexEnvValue[n] < 0;

            if (growthStopped || outsideGrowthField)
            {
                isVertexFixed[n]         = true;
                perVertexGrowthRate[n]   = 0;
                IGoal anchor             = new AnchorXYZ(updatedPositions[n], true, true, true, 100);
                anchor.PIndex            = new int[1] { n };
                anchorGoals.Add(anchor);
            }

            halfedgeMesh.Vertices.SetVertex(n, updatedPositions[n]);
        }
    }

    // Outputs 
    A = halfedgeMesh;
    B = perVertexGrowthRate;
    C = _iterationCount++;
    E = perVertexEnvValue;
    F = isVertexFixed;
}


// Persistent State Fields

List<double> perVertexTargetLength;
List<double> perVertexGrowthRate;
List<double> perVertexEnvValue = new List<double>();
List<bool>   isVertexFixed;
PlanktonMesh halfedgeMesh;
Mesh         collideMesh;
int          _iterationCount = 0;

// Environment data grid
double[][] envData;
int        xCount, yCount;
UVInterval uvDomain;

// Physics
KangarooSolver.PhysicalSystem physicsSystem;
List<IGoal> springGoals;
List<IGoal> hingeGoals;
List<IGoal> anchorGoals;


// Bilinear Interpolation of Environment Data 
public double Interpolate(Point3d point)
{
    double xNorm = xCount * (point.X - uvDomain.U0) / uvDomain.U.Length;
    int    xIdx  = xNorm > 0 ? (xNorm < xCount - 2 ? (int) xNorm : xCount - 2) : 0;

    double yNorm = yCount * (point.Y - uvDomain.V0) / uvDomain.V.Length;
    int    yIdx  = yNorm > 0 ? (yNorm < yCount - 2 ? (int) yNorm : yCount - 2) : 0;

    double xFrac = xNorm - xIdx; // fractional part in X
    double yFrac = yNorm - yIdx; // fractional part in Y

    double bottomRow = envData[xIdx][yIdx]     * (1 - xFrac) + envData[xIdx + 1][yIdx]     * xFrac;
    double topRow    = envData[xIdx][yIdx + 1] * (1 - xFrac) + envData[xIdx + 1][yIdx + 1] * xFrac;

    return bottomRow * (1 - yFrac) + topRow * yFrac;
}


// Hinge Vertex Lookup
public int[] GetHingeVertices(int halfedgeIndex, PlanktonMesh mesh, out double restAngle)
{
    int[] verts = new int[4];

    // v0: start of the halfedge
    verts[0] = mesh.Halfedges[halfedgeIndex].StartVertex;
    Point3d P0 = mesh.Vertices[verts[0]].ToPoint3d();

    // v1: start of the pair (other end of edge)
    int pairHe = mesh.Halfedges.GetPairHalfedge(halfedgeIndex);
    verts[1] = mesh.Halfedges[pairHe].StartVertex;
    Point3d P1 = mesh.Vertices[verts[1]].ToPoint3d();

    // v2: tip of the triangle on the pair side
    int pairPrevHe = mesh.Halfedges[pairHe].PrevHalfedge;
    verts[2] = mesh.Halfedges[pairPrevHe].StartVertex;
    Point3d P2 = mesh.Vertices[verts[2]].ToPoint3d();

    // v3: tip of the triangle on the primary side
    int primaryPrevHe = mesh.Halfedges[halfedgeIndex].PrevHalfedge;
    verts[3] = mesh.Halfedges[primaryPrevHe].StartVertex;
    Point3d P3 = mesh.Vertices[verts[3]].ToPoint3d();

    // Dihedral angle between the two adjacent triangles
    Vector3d edgeVec    = P1 - P0;
    Vector3d toV2       = P2 - P0;
    Vector3d toV3       = P3 - P0;
    Vector3d normalLeft  = Vector3d.CrossProduct(toV2, edgeVec);
    Vector3d normalRight = Vector3d.CrossProduct(edgeVec, toV3);

    double angle = Vector3d.VectorAngle(normalLeft, normalRight, new Plane(P0, edgeVec));
    if (angle > Math.PI) angle -= 2 * Math.PI;
    restAngle = angle;

    return verts;
}


// GrowthPush Goal
// Pushes vertices apart when they come within a threshold distance of each other.
// Uses a tighter distance for mesh-adjacent pairs and a looser one for non-adjacent pairs.

public class GrowthPush : GoalObject
{
    public double Strength = 1.0;
    public double sqrPushDistAdjacent;
    public double sqrPushDistNonAdjacent;
    public double pushDistAdjacent;
    public double pushDistNonAdjacent;
    PlanktonMesh sourceMesh;

    private bool AreVerticesAdjacent(int vertA, int vertB)
    {
        IEnumerable<int> neighborHalfedges = sourceMesh.Halfedges.GetVertexCirculator(
            sourceMesh.Vertices[vertA].OutgoingHalfedge);

        foreach (var he in neighborHalfedges)
        {
            if (sourceMesh.Halfedges[sourceMesh.Halfedges[he].NextHalfedge].StartVertex == vertB)
                return true;
        }
        return false;
    }

    public GrowthPush(double distAdjacent, double distNonAdjacent, PlanktonMesh mesh)
    {
        PIndex   = Enumerable.Range(0, mesh.Vertices.Count).ToArray();
        int num  = PIndex.Length;
        Move     = new Vector3d[num];
        Weighting = new double[num];

        for (int i = 0; i < num; i++)
            Weighting[i] = Strength;

        sqrPushDistAdjacent    = distAdjacent * distAdjacent;
        sqrPushDistNonAdjacent = distNonAdjacent * distNonAdjacent;
        pushDistAdjacent       = distAdjacent;
        pushDistNonAdjacent    = distNonAdjacent;
        sourceMesh             = mesh;
    }

    public override void Calculate(List<KangarooSolver.Particle> particles)
    {
        int num = PIndex.Length;

        // Sort vertices by Z so we can early-exit inner loop when Z gap exceeds push distance
        double[] zCoords = new double[num];
        for (int i = 0; i < num; i++)
            zCoords[i] = particles[PIndex[i]].Position.Z;

        Array.Sort(zCoords, PIndex);

        Parallel.For(0, PIndex.Length - 1, i =>
        {
            for (int j = 1; (i + j) < PIndex.Length; j++)
            {
                int k     = i + j;
                int idxHi = PIndex[k];
                int idxLo = PIndex[i];

                bool adjacent   = AreVerticesAdjacent(idxHi, idxLo);
                double pushDist = adjacent ? pushDistAdjacent : pushDistNonAdjacent;
                double sqrDist  = adjacent ? sqrPushDistAdjacent : sqrPushDistNonAdjacent;

                Vector3d separation = particles[idxHi].Position - particles[idxLo].Position;

                if (separation.Z >= pushDist) break; // further vertices in sorted order are too far

                double sqrLen = separation.SquareLength;
                if (sqrLen < sqrDist)
                {
                    double len    = Math.Sqrt(sqrLen);
                    double factor = 1.0 - pushDist / len;
                    Vector3d push = 0.1 * separation * factor;
                    Move[i] += push;
                    Move[k] -= push;
                }
            }
        });
    }
}


// ── MeshPush Goal ─────────────────────────────────────────────────────────────
// Pushes vertices away from a collision mesh surface when within a threshold distance.

public class MeshPush : GoalObject
{
    public double Strength = 1.0;
    PlanktonMesh sourceMesh;
    Mesh         targetMesh;
    double       pushDist;
    double       sqrPushDist;

    public MeshPush(PlanktonMesh meshParticles, Mesh collisionTarget, double dist)
    {
        PIndex    = Enumerable.Range(0, meshParticles.Vertices.Count).ToArray();
        int num   = PIndex.Length;
        Move      = new Vector3d[num];
        Weighting = new double[num];

        for (int i = 0; i < num; i++)
            Weighting[i] = Strength;

        sourceMesh  = meshParticles;
        targetMesh  = collisionTarget;
        sqrPushDist = dist * dist;
        pushDist    = dist;
    }

    public override void Calculate(List<KangarooSolver.Particle> particles)
    {
        Parallel.For(0, PIndex.Length, i =>
        {
            Point3d vertexPos      = particles[PIndex[i]].Position;
            var     closestPoint   = targetMesh.ClosestMeshPoint(vertexPos, 1000);
            Vector3d toSurface     = closestPoint.Point - vertexPos;
            double   sqrLen        = toSurface.SquareLength;

            if (sqrLen < sqrPushDist)
            {
                double   len    = Math.Sqrt(sqrLen);
                double   factor = 1.0 - pushDist / len;
                Move[i] = 0.5 * toSurface * factor;
            }
        });
    }
  // </Custom additional code> 
}
