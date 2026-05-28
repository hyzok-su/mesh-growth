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
  private void RunScript(int xnum, int ynum, List<double> data, UVInterval interval2d, Mesh pmesh, Mesh colmesh, double splitratio, double decsentrate, List<double> targetlengths, List<double> growthrate, double growthdecsent, List<bool> isanchored, double hingestrength, bool reset, int viewnum, ref object A, ref object B, ref object C, ref object D, ref object E, ref object F)
  {

    //pre
    bool stop;
    double hingestr = hingestrength;
    double springstr = 1;
    double dist1 = 1;
    double dist2 = 1.5;
    double sr = splitratio >= 0 ? (splitratio <= 1 ? splitratio : 1) : 0;


    if(reset)
    {
      //edata
      edata = new double[xnum][];
      for(int i = 0;i < xnum;i++)
      {
        edata[i] = new double[ynum];
        for(int j = 0;j < ynum;j++)
        {
          edata[i][j] = data[j + i * ynum];
        }
      }
      xNum = xnum;yNum = ynum;
      uvDom = interval2d;

      colMesh = colmesh;

      tL = targetlengths;
      P = pmesh.ToPlanktonMesh();
      gR = growthrate;

      PSys = new KangarooSolver.PhysicalSystem();
      anchors = new List<IGoal>();
      springs = new List<IGoal>();
      hinges = new List<IGoal>();
      skip = new List<bool>();
      _counter = 0;
      int cc = 0;
      foreach(var v in P.Vertices)
      {
        var pt = v.ToPoint3d();
        var ee = Interpolate(pt);
        eR.Add(ee);
        PSys.AddParticle(pt, 10);
        if (!isanchored[cc] && ee > 0)
        {
          skip.Add(false);
        }
        else
        {
          IGoal anch = new AnchorXYZ(pt, true, true, true, 100);
          anch.PIndex = new int[1]{cc};
          anchors.Add(anch);
          skip.Add(true);
        }
        cc++;
      }

      for(int i = 0;i < P.Halfedges.Count / 2;i++)
      {
        var h = P.Halfedges[2 * i];
        var ph = P.Halfedges[2 * i + 1];
        int s = h.StartVertex;
        int e = ph.StartVertex;
        double length_now = (P.Vertices[s].ToPoint3d() - P.Vertices[e].ToPoint3d()).Length;
        IGoal spr = new Spring(s, e, length_now, springstr);
        springs.Add(spr);
        double a;
        int[] pts = GetHinge(2 * i, P, out a);
        if(h.AdjacentFace != -1 && ph.AdjacentFace != -1)
        {
          IGoal hin = new Hinge(pts[0], pts[1], pts[2], pts[3], 0, hingestr);
          hinges.Add(hin);
        }
        else
        {
          IGoal hin = new Hinge(pts[0], pts[1], pts[2], pts[3], 0, 0);
          hinges.Add(hin);
        }
      }
    }

    double t = 1;
    var counter = 0;
    var dr = decsentrate; dr = dr < 0.96 ? (dr < 0 ? 0 : dr) : 0.96;
    var gd = growthdecsent; gd = gd < 1 ? (gd < 0.6 ? 0.6 : gd) : 1;

    //run
    if(!reset)
    {
      //remesh start
      do
      {
        stop = true;
        var hs = P.Halfedges;
        var ec = hs.Count / 2;
        var el = hs.GetLengths();
        var vs = P.Vertices;
        for(int i = 0;i < ec;i++)
        {
          if (hs[2 * i].IsUnused == false)
          {
            int vStart = hs[2 * i].StartVertex;
            int vEnd = hs[2 * i + 1].StartVertex;
            //only when skip is false it executes splitting or collapse.
            //there are two reasons that it skips.
            if (!skip[vStart] || !skip[vEnd])
            {
              var tl = (tL[vStart] + tL[vEnd]) / 2;
              //split too long
              if(el[2 * i] > (1 + t) * 1.25f * tl)
              {
                int spe = hs.TriangleSplitEdge(2 * i);
                if(spe != -1)//succeed
                {
                  var g1 = gR[vStart] *= gd;
                  var g2 = gR[vEnd] *= gd;
                  var g3 = g1 > g2 ? g2 : g1;//minimum
                  var g4 = (g1 + g2) / 2;//average
                  var g5 = (1 - sr) * g3 + sr * g4;//somewhere between minimum and average
                  stop = false;
                  int spv = hs[spe].StartVertex;
                  skip.Add(false);
                  tL.Add(tl);
                  gR.Add(g5);

                  //add particle
                  var spt = vs[vStart].ToPoint3d();
                  var ept = vs[vEnd].ToPoint3d();
                  var mid = (spt + ept) / 2;
                  PSys.AddParticle(mid, 10);
                  eR.Add(Interpolate(mid));

                  //resize goals
                  int sizediff = hs.Count / 2 - springs.Count;
                  springs.AddRange(new Spring[sizediff]);
                  hinges.AddRange(new Hinge[sizediff]);
                  //loop add new & mark skip
                  IEnumerable <int> nhs = hs.GetVertexCirculator(spe);D = hs[spe].StartVertex;
                  foreach(var nh in nhs)
                  {
                    int edge = nh / 2;
                    double restangle;
                    int[] hpts = GetHinge(2 * edge, P, out restangle);
                    //hinge
                    if(P.Halfedges[nh].AdjacentFace != -1 && P.Halfedges[P.Halfedges.GetPairHalfedge(nh)].AdjacentFace != -1)
                    {
                      IGoal hin = new Hinge(hpts[0], hpts[1], hpts[2], hpts[3], 0, hingestr);
                      hinges[edge] = hin;
                    }
                    else
                    {
                      IGoal hin = new Hinge(hpts[0], hpts[1], hpts[2], hpts[3], 0, 0);
                      hinges[edge] = hin;
                    }

                    //spring
                    var lengthhere = vs[hpts[0]].ToPoint3d().DistanceTo(vs[hpts[1]].ToPoint3d());
                    IGoal spr = new Spring(hpts[0], hpts[1], lengthhere, springstr);
                    springs[edge] = spr;
                  }
                }
              }

            }
          }
        }
        P.Compact();

        hs = P.Halfedges;
        ec = hs.Count / 2;
        el = hs.GetLengths();
        vs = P.Vertices;
        //merge too short
        for(int i = 0;i < ec;i++)
        {
          if (hs[2 * i].IsUnused == false)
          {
            int vStart = hs[2 * i].StartVertex;
            int vEnd = hs[2 * i + 1].StartVertex;
            if (!skip[vStart] || !skip[vEnd])
            {
              var tl = (tL[vStart] + tL[vEnd]) / 2;
              if(el[2 * i] < (1 - t) * 0.25f * tl)
              {
                //store the index of edge to remove
                int next = hs[2 * i].NextHalfedge;
                int prev = hs[2 * i + 1].PrevHalfedge;
                int left = hs[2 * i].AdjacentFace;
                int right = hs[2 * i + 1].AdjacentFace;

                int clStart = hs[2 * i].StartVertex;
                int clEnd = hs[2 * i + 1].StartVertex;
                float x = (vs[clStart].X + vs[clEnd].X) / 2;
                float y = (vs[clStart].Y + vs[clEnd].Y) / 2;
                float z = (vs[clStart].Z + vs[clEnd].Z) / 2;
                var g1 = gR[clStart];
                var g2 = gR[clEnd];
                var g = g1 > g2 ? g1 : g2;
                int cle = hs.CollapseEdge(2 * i);
                if(cle != -1)
                {
                  stop = false;
                  vs.SetVertex(clStart, x, y, z);
                  tL[clStart] = tl;
                  gR[clStart] = g;
                  tL.RemoveAt(clEnd);
                  gR.RemoveAt(clEnd);

                  //remove
                  springs.RemoveAt(i);
                  hinges.RemoveAt(i);
                  PSys.DeleteParticle(clEnd);
                  if(left != -1)
                  {
                    springs.RemoveAt(next / 2);
                    hinges.RemoveAt(next / 2);
                  }
                  if(right != -1)
                  {
                    springs.RemoveAt(prev / 2);
                    hinges.RemoveAt(prev / 2);
                  }
                  //loop change p & mark skip
                  IEnumerable <int> nhs = hs.GetVertexCirculator(cle);
                  foreach(var nh in nhs)
                  {
                    int edge = nh / 2;
                    //hinge
                    hinges[edge].PIndex[0] = hinges[edge].PIndex[0] == clEnd ? clStart : hinges[edge].PIndex[0];
                    hinges[edge].PIndex[1] = hinges[edge].PIndex[1] == clEnd ? clStart : hinges[edge].PIndex[1];
                    hinges[edge].PIndex[2] = hinges[edge].PIndex[2] == clEnd ? clStart : hinges[edge].PIndex[2];
                    hinges[edge].PIndex[3] = hinges[edge].PIndex[3] == clEnd ? clStart : hinges[edge].PIndex[3];
                    //spring
                    springs[edge].PIndex[0] = springs[edge].PIndex[0] == clEnd ? clStart : springs[edge].PIndex[0];
                    springs[edge].PIndex[1] = springs[edge].PIndex[1] == clEnd ? clStart : springs[edge].PIndex[1];
                    //skip
                    skip[hs[hs.GetPairHalfedge(nh)].StartVertex] = true;
                  }
                }
              }
            }
          }
        }
        P.Compact();
        ec = hs.Count / 2;
        for(int i = 0; i < ec;i++)
        {
          if(!hs[2 * i].IsUnused
            && hs[2 * i].AdjacentFace != -1
            && hs[2 * i + 1].AdjacentFace != -1
            && hs[hs.GetPairHalfedge(hs[2 * i].NextHalfedge)].AdjacentFace != -1
            && hs[hs.GetPairHalfedge(hs[2 * i].PrevHalfedge)].AdjacentFace != -1
            && hs[hs.GetPairHalfedge(hs[2 * i + 1].NextHalfedge)].AdjacentFace != -1
            && hs[hs.GetPairHalfedge(hs[2 * i + 1].PrevHalfedge)].AdjacentFace != -1
            )
          {
            int v1 = hs[2 * i].StartVertex;
            int v2 = hs[2 * i + 1].StartVertex;
            int v3 = hs[hs[hs[2 * i].NextHalfedge].NextHalfedge].StartVertex;
            int v4 = hs[hs[hs[2 * i + 1].NextHalfedge].NextHalfedge].StartVertex;
            Point3d P1 = vs[v1].ToPoint3d();
            Point3d P2 = vs[v2].ToPoint3d();
            Point3d P3 = vs[v3].ToPoint3d();
            Point3d P4 = vs[v4].ToPoint3d();
            double A1 = Vector3d.VectorAngle(new Vector3d(P3 - P1), new Vector3d(P4 - P1))
              + Vector3d.VectorAngle(new Vector3d(P4 - P2), new Vector3d(P3 - P2));
            double A2 = Vector3d.VectorAngle(new Vector3d(P1 - P4), new Vector3d(P2 - P4))
              + Vector3d.VectorAngle(new Vector3d(P2 - P3), new Vector3d(P1 - P3));
            if (A2 > A1)
            {
              stop = false;
              P.Halfedges.FlipEdge(2 * i);

              //modify springs
              springs[i] = new Spring(v3, v4, P3.DistanceTo(P4), springstr);
              //modify hinges
              double a;
              GetHinge(2 * i, P, out a);
              hinges[i] = new Hinge(v3, v4, v2, v1, 0, hingestr);
            }
          }
        }
        stop = t > 0.05 ? false : stop;
        t *= dr;
        counter++;
      }while(!(stop) && counter < 101);Print("{0}", counter - 1);
      P.Compact();

      //growth start
      var goals = new List<IGoal>();
      for(int i = 0;i < P.Halfedges.Count / 2;i++)
      {
        Spring s = springs[i] as Spring;
        s.RestLength += (gR[s.PIndex[0]] + gR[s.PIndex[1]] ) * (eR[s.PIndex[0]] + eR[s.PIndex[1]] ) / 4; goals.Add((IGoal) s);
        Hinge hge = hinges[i] as Hinge;
        if(hge.Strength != 0){goals.Add((IGoal) hge);}
      }
      goals.AddRange(anchors);
      goals.Add(new GrowthPush(dist1, dist2, P));
      goals.Add(new MeshPush(P, colMesh, 0.5));
      PSys.Step(goals, true, 0.1);
      var particles = PSys.GetPositionArray();
      for(int n = 0;n < particles.Count();n++)
      {
        Vector3d sp = particles[n] - P.Vertices[n].ToPoint3d();
        eR[n] = Interpolate(particles[n]);
        //when growthrate<0.0001, stop growth
        if(!skip[n] && sp.Length < 0.0001 && gR[n] < 0.0001 || eR[n] < 0)
        {
          skip[n] = true;
          gR[n] = 0;
          IGoal anch = new AnchorXYZ(particles[n], true, true, true, 100);
          anch.PIndex = new int[1]{n};
          anchors.Add(anch);
        }
        P.Vertices.SetVertex(n, particles[n]);
      }
    }
    A = P;
    B = gR;
    C = _counter++;
    E = eR;
    F = skip;


  }

  // <Custom additional code> 



  List<double> tL;//target length
  List<double> gR;//growth rate
  List<double> eR = new List<double>();
  List<bool> skip;
  PlanktonMesh P;
  Mesh colMesh;
  int _counter = 0;

  //edata
  double[][] edata;
  int xNum,yNum;
  UVInterval uvDom;
  //interpolaton
  public double Interpolate(Point3d pt)
  {
    var x = xNum * (pt.X - uvDom.U0) / uvDom.U.Length;
    int xx = x > 0 ? (x < xNum - 2 ? (int) x : xNum - 2) : 0;

    var y = yNum * (pt.Y - uvDom.V0) / uvDom.V.Length;
    int yy = y > 0 ? (y < yNum - 2 ? (int) y : yNum - 2) : 0;


    var y1 =
      edata[xx][yy] * (1 - x + xx)
      +
      edata[xx + 1][yy] * (x - xx);

    var y2 =
      edata[xx][yy + 1] * (1 - x + xx)
      +
      edata[xx + 1][yy + 1] * (x - xx);

    return
      y1 * (1 - y + yy)
      +
      y2 * (y - yy);
  }

  //physics
  KangarooSolver.PhysicalSystem PSys;
  List<IGoal> springs;
  List<IGoal> hinges;
  List<IGoal> anchors;

  public int[] GetHinge(int h, PlanktonMesh p, out double a)
  {
    int[] output = new int[4];

    output[0] = p.Halfedges[h].StartVertex;
    Point3d P0 = p.Vertices[p.Halfedges[h].StartVertex].ToPoint3d();

    int ph = p.Halfedges.GetPairHalfedge(h);
    output[1] = p.Halfedges[ph].StartVertex;
    Point3d P1 = p.Vertices[p.Halfedges[ph].StartVertex].ToPoint3d();

    int php = p.Halfedges[ph].PrevHalfedge;
    output[2] = p.Halfedges[php].StartVertex;
    Point3d P2 = p.Vertices[p.Halfedges[php].StartVertex].ToPoint3d();

    int hp = p.Halfedges[h].PrevHalfedge;
    output[3] = p.Halfedges[hp].StartVertex;
    Point3d P3 = p.Vertices[p.Halfedges[hp].StartVertex].ToPoint3d();

    Vector3d V01 = P1 - P0;
    Vector3d V02 = P2 - P0;
    Vector3d V03 = P3 - P0;

    Vector3d Cross0 = Vector3d.CrossProduct(V02, V01);
    Vector3d Cross1 = Vector3d.CrossProduct(V01, V03);

    double temp = Vector3d.VectorAngle(Cross0, Cross1, new Plane(P0, V01));
    if (temp > Math.PI) { temp -= 2 * Math.PI; }
    a = temp;
    return output;
  }




  //this class of instance refresh in every iteration.
  public class GrowthPush:GoalObject
  {
    public double Strength = 1;
    public double SqrPushDist1;
    public double SqrPushDist2;
    public double PushDist1;
    public double PushDist2;
    PlanktonMesh Pm;


    private bool CheckNeighbor(int v1, int v2)
    {
      IEnumerable <int> nhs = Pm.Halfedges.GetVertexCirculator(Pm.Vertices[v1].OutgoingHalfedge);
      foreach( var nh in nhs)
      {
        if (Pm.Halfedges[Pm.Halfedges[nh].NextHalfedge].StartVertex == v2)
        {return true;}
      }
      return false;
    }


    //GrowthPush
    public GrowthPush(double dist1, double dist2, PlanktonMesh pm)
    {
      //base
      PIndex = Enumerable.Range(0, pm.Vertices.Count).ToArray();
      int num = PIndex.Length;
      Move = new Vector3d[num];
      Weighting = new double[num];
      for (int i = 0; i < num; i++)
      {
        Weighting[i] = Strength;
      }
      //this
      SqrPushDist1 = dist1 * dist1;
      SqrPushDist2 = dist2 * dist2;
      PushDist1 = dist1;
      PushDist2 = dist2;
      Pm = pm;
    }

    public override void Calculate(List < KangarooSolver.Particle > p)
    {
      int num = PIndex.Length;
      double[] Zcoord = new double[num];

      for (int i = 0; i < num; i++)
      {
        Zcoord[i] = p[PIndex[i]].Position.Z;
        Move[i] = Vector3d.Zero;
      }
      Array.Sort(Zcoord, PIndex);


      Parallel.For(0, PIndex.Length - 1, i =>
        {
        for (int j = 1; (i + j) < PIndex.Length; j++)
        {
          int k = i + j;
          int ind1 = PIndex[k];
          int ind2 = PIndex[i];
          bool isPair = CheckNeighbor(ind1, ind2);
          Vector3d Sp = p[ind1].Position - p[ind2].Position;
          if(isPair)
          {
            if (Sp.Z < PushDist1)
            {
              for(var sl = Sp.SquareLength;sl < SqrPushDist1;)
              {
                var l = Math.Sqrt(sl);
                double factor = 1 - PushDist1 / l;
                Vector3d dmove = 0.1 * Sp * factor;
                Move[i] += dmove;
                Move[k] -= dmove;
                break;
              }
            }
            else { break; }
          }
          else
          {
            if (Sp.Z < PushDist2)
            {
              for(var sl = Sp.SquareLength;sl < SqrPushDist2;)
              {
                var l = Math.Sqrt(sl);
                double factor = 1 - PushDist2 / l;
                Vector3d dmove = 0.1 * Sp * factor;
                Move[i] += dmove;
                Move[k] -= dmove;
                break;
              }
            }
            else { break; }
          }
        }
        });

    }
  }
  //MeshPush
  public class MeshPush : GoalObject
  {

    public double Strength = 1;
    PlanktonMesh Pm;
    Mesh M;
    double Dist;
    double SquareDist;

    public MeshPush(PlanktonMesh pm, Mesh m, double dist)
    {
      //base
      PIndex = Enumerable.Range(0, pm.Vertices.Count).ToArray();
      int num = PIndex.Length;
      Move = new Vector3d[num];
      Weighting = new double[num];
      for (int i = 0; i < num; i++)
      {
        Weighting[i] = Strength;
      }
      //this
      Pm = pm;M = m;SquareDist = dist * dist;Dist = dist;
    }

    public override void Calculate(List < KangarooSolver.Particle > p)
    {
      Parallel.For(0, PIndex.Length, i =>
        {
        Point3d ThisPt = p[PIndex[i]].Position;
        var MP = M.ClosestMeshPoint(ThisPt, 1000);
        var Push = MP.Point - ThisPt;
        var sl = Push.SquareLength;
        if(sl < SquareDist)
        {
          var l = Math.Sqrt(sl);
          double factor = 1 - Dist / l;
          Vector3d dmove = 0.5 * Push * factor;
          Move[i] = dmove;
        }
        });
    }
  }

  // </Custom additional code> 
}
