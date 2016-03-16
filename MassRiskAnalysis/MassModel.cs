using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Text;

namespace MassRiskAnalysis
{
    public class MassNodeList : List<MassNode>
    {
        public bool IsEmpty
        {
            get { return this.Count == 0 ? true : false; }
        }
        public static void Update(List<MassNode> node_list)
        {
            Queue<MassNode> temp_list = new Queue<MassNode>();
            Parallel.ForEach(node_list, (node) =>
            {
                node.level = node.ChildNodes.Count;
                if (node.level == 0)
                {
                    temp_list.Enqueue(node);
                    if (!node.Updated && node.ParentNode != null)
                    {
                        node.ParentNode.Updated = false;   
                    }
                    node.Updated = true;
                }
            });
            while (temp_list.Count > 0)
            {
                MassNode node = temp_list.Dequeue();
                if (node.ParentNode != null)
                {
                    node.ParentNode.level--;
                    if (node.ParentNode.level == 0)
                        temp_list.Enqueue(node.ParentNode);
                }
                if (node.Updated) continue;
                double M = 0;
                double MX = 0, MY = 0, MZ = 0;
                double MXX = 0, MYY = 0, MZZ = 0, MXY = 0, MYZ = 0, MZX = 0;
                double IXX = 0, IYY = 0, IZZ = 0, IXY = 0, IYZ = 0, IZX = 0;
                Parallel.ForEach(node.ChildNodes, (cnode) =>
                {
                    double Mi = cnode.Mass;
                    double MXi = Mi * cnode.Xcg;
                    double MYi = Mi * cnode.Ycg;
                    double MZi = Mi * cnode.Zcg;
                    double MXXi = MXi * cnode.Xcg;
                    double MYYi = MYi * cnode.Ycg;
                    double MZZi = MZi * cnode.Zcg;
                    double MXYi = MXi * cnode.Ycg;
                    double MYZi = MYi * cnode.Zcg;
                    double MZXi = MZi * cnode.Xcg;
                    M += Mi;
                    MX += MXi;
                    MY += MYi;
                    MZ += MZi;
                    MXX += MXXi;
                    MYY += MYYi;
                    MZZ += MZZi;
                    MXY += MXYi;
                    MYZ += MYZi;
                    MZX += MZXi;
                    IXX += cnode.Ixx;
                    IYY += cnode.Iyy;
                    IZZ += cnode.Izz;
                    IXY += cnode.Ixy;
                    IYZ += cnode.Iyz;
                    IZX += cnode.Izx;
                });
                if (Math.Abs(M) <= double.Epsilon)
                {
                    node.ResetValue();
                }
                else
                {
                    double mass = M;
                    double xcg = MX / M;
                    double ycg = MY / M;
                    double zcg = MZ / M;
                    double xxcg = xcg * xcg;
                    double yycg = ycg * ycg;
                    double zzcg = zcg * zcg;
                    double xycg = xcg * ycg;
                    double yzcg = ycg * zcg;
                    double zxcg = zcg * xcg;
                    double ixx = IXX + MYY + MZZ + (yycg + zzcg) * M - 2 * ycg * MY - 2 * zcg * MZ;
                    double iyy = IYY + MZZ + MXX + (zzcg + xxcg) * M - 2 * zcg * MZ - 2 * xcg * MX;
                    double izz = IZZ + MXX + MYY + (xxcg + yycg) * M - 2 * xcg * MX - 2 * ycg * MY;
                    double ixy = IXY + xycg * M - xcg * MY - ycg * MX + MXY;
                    double iyz = IYZ + yzcg * M - ycg * MZ - zcg * MY + MYZ;
                    double izx = IZX + zxcg * M - zcg * MX - xcg * MZ + MZX;
                    node.SetValue(mass, xcg, ycg, zcg, ixx, iyy, izz, ixy, iyz, izx);
                }
                node.Updated = true;
            }
        }
    }

    public class MassNodeListForSample : List<MassNodeForSample>
    {
        public static void GenerateSample(List<MassNodeForSample> node_list, int n)
        {
            Queue<MassNodeForSample> temp_list = new Queue<MassNodeForSample>();
            foreach (MassNodeForSample node in node_list)
            {
                if (node.ChildNodes.IsEmpty)
                {
                    temp_list.Enqueue(node);
                    node.SampleList.Clear();
                    node.SampleList = new List<MassProperties>(n);
                    Parallel.ForEach(node.SampleList, (sample) =>
                    {
                        sample = new MassProperties(node.SampleModel.Sample(),
                            node.Xcg, node.Ycg, node.Zcg,
                            node.Ixx, node.Iyy, node.Izz,
                            node.Ixy, node.Iyz, node.Izx);
                    });
                }
            };
            while (temp_list.Count > 0)
            {
                MassNodeForSample node = temp_list.Dequeue();
                if (node.ParentNode != null)
                {
                    node.ParentNode.level--;
                    if (node.ParentNode.level == 0)
                        temp_list.Enqueue((MassNodeForSample)node.ParentNode);
                }
                if (node.ChildNodes.IsEmpty) return;
                node.SampleList.Clear();
                node.SampleList = new List<MassProperties>(n);
                Parallel.ForEach(node.SampleList, (sample) =>
                {
                    sample = new MassProperties(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                    int i = node.SampleList.IndexOf(sample);
                    double M = 0;
                    double MX = 0, MY = 0, MZ = 0;
                    double MXX = 0, MYY = 0, MZZ = 0, MXY = 0, MYZ = 0, MZX = 0;
                    double IXX = 0, IYY = 0, IZZ = 0, IXY = 0, IYZ = 0, IZX = 0;
                    foreach(MassNodeForSample snode in node.ChildNodes)
                    {
                        double Mi = snode.SampleList[i].Mass;
                        double MXi = Mi * snode.SampleList[i].Xcg;
                        double MYi = Mi * snode.SampleList[i].Ycg;
                        double MZi = Mi * snode.SampleList[i].Zcg;
                        double MXXi = MXi * snode.SampleList[i].Xcg;
                        double MYYi = MYi * snode.SampleList[i].Ycg;
                        double MZZi = MZi * snode.SampleList[i].Zcg;
                        double MXYi = MXi * snode.SampleList[i].Ycg;
                        double MYZi = MYi * snode.SampleList[i].Zcg;
                        double MZXi = MZi * snode.SampleList[i].Xcg;
                        M += Mi;
                        MX += MXi;
                        MY += MYi;
                        MZ += MZi;
                        MXX += MXXi;
                        MYY += MYYi;
                        MZZ += MZZi;
                        MXY += MXYi;
                        MYZ += MYZi;
                        MZX += MZXi;
                        IXX += snode.SampleList[i].Ixx;
                        IYY += snode.SampleList[i].Iyy;
                        IZZ += snode.SampleList[i].Izz;
                        IXY += snode.SampleList[i].Ixy;
                        IYZ += snode.SampleList[i].Iyz;
                        IZX += snode.SampleList[i].Izx;
                    };
                    if (Math.Abs(M) <= double.Epsilon)
                    {
                        sample.ResetValue();
                    }
                    else
                    {
                        double mass = M;
                        double xcg = MX / M;
                        double ycg = MY / M;
                        double zcg = MZ / M;
                        double xxcg = xcg * xcg;
                        double yycg = ycg * ycg;
                        double zzcg = zcg * zcg;
                        double xycg = xcg * ycg;
                        double yzcg = ycg * zcg;
                        double zxcg = zcg * xcg;
                        double ixx = IXX + MYY + MZZ + (yycg + zzcg) * M - 2 * ycg * MY - 2 * zcg * MZ;
                        double iyy = IYY + MZZ + MXX + (zzcg + xxcg) * M - 2 * zcg * MZ - 2 * xcg * MX;
                        double izz = IZZ + MXX + MYY + (xxcg + yycg) * M - 2 * xcg * MX - 2 * ycg * MY;
                        double ixy = IXY + xycg * M - xcg * MY - ycg * MX + MXY;
                        double iyz = IYZ + yzcg * M - ycg * MZ - zcg * MY + MYZ;
                        double izx = IZX + zxcg * M - zcg * MX - xcg * MZ + MZX;
                        sample.SetValue(mass, xcg, ycg, zcg, ixx, iyy, izz, ixy, iyz, izx);
                    }
                });

            }

        }
    }

    public class MassNode
    {
        private bool _updated = false;
        public bool Updated { get { return _updated; } set { _updated = value; } }
        private MassProperties Properties = new MassProperties(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        public double Mass { get { return this.Properties.Mass; } }
        public double Xcg { get { return this.Properties.Xcg; } }
        public double Ycg { get { return this.Properties.Ycg; } }
        public double Zcg { get { return this.Properties.Zcg; } }
        public double Ixx { get { return this.Properties.Ixx; } }
        public double Iyy { get { return this.Properties.Iyy; } }
        public double Izz { get { return this.Properties.Izz; } }
        public double Ixy { get { return this.Properties.Ixy; } }
        public double Iyz { get { return this.Properties.Iyz; } }
        public double Izx { get { return this.Properties.Izx; } }
        public int level = 0;
        public MassNode ParentNode = null;

        public MassNodeList ChildNodes = new MassNodeList();

        public void SetValue(double mass = 0, double xcg = 0, double ycg = 0, double zcg = 0,
            double ixx = 0, double iyy = 0, double izz = 0, double ixy = 0, double iyz = 0, double izx = 0)
        {
            this.Properties.SetValue(mass, xcg, ycg, zcg, ixx, iyy, izz, ixy, iyz, izx);
        }

        public void SetValueNew(double mass = 0, double xcg = 0, double ycg = 0, double zcg = 0,
    double ixx = 0, double iyy = 0, double izz = 0, double ixy = 0, double iyz = 0, double izx = 0)
        {
            this.Properties.SetValue(mass, xcg, ycg, zcg, ixx, iyy, izz, ixy, iyz, izx);
            MassNode node = this;
            while (node.Updated)
            {
                node.Updated = false;
                if (node.ParentNode != null)
                    node = node.ParentNode;
            }
        }


        public void ResetValue()
        {
            this.Properties.ResetValue();
        }

        public void ResetValueNew()
        {
            this.Properties.ResetValue();
            MassNode node = this;
            while (node.Updated)
            {
                node.Updated = false;
                if (node.ParentNode != null)
                    node = node.ParentNode;
            }
        }

        //public void Update()
        //{
        //    if (this.ChildNodes.IsEmpty)
        //    {
        //        this._updated = true;
        //        return;
        //    }
        //    double M = 0;
        //    double MX = 0, MY = 0, MZ = 0;
        //    double MXX = 0, MYY = 0, MZZ = 0, MXY = 0, MYZ = 0, MZX = 0;
        //    double IXX = 0, IYY = 0, IZZ = 0, IXY = 0, IYZ = 0, IZX = 0;
        //    Parallel.ForEach(ChildNodes, (node) =>
        //    {
        //        node.Update();
        //        double Mi = node.Properties.Mass;
        //        double MXi = Mi * node.Properties.Xcg;
        //        double MYi = Mi * node.Properties.Ycg;
        //        double MZi = Mi * node.Properties.Zcg;
        //        double MXXi = MXi * node.Properties.Xcg;
        //        double MYYi = MYi * node.Properties.Ycg;
        //        double MZZi = MZi * node.Properties.Zcg;
        //        double MXYi = MXi * node.Properties.Ycg;
        //        double MYZi = MYi * node.Properties.Zcg;
        //        double MZXi = MZi * node.Properties.Xcg;
        //        M += Mi;
        //        MX += MXi;
        //        MY += MYi;
        //        MZ += MZi;
        //        MXX += MXXi;
        //        MYY += MYYi;
        //        MZZ += MZZi;
        //        MXY += MXYi;
        //        MYZ += MYZi;
        //        MZX += MZXi;
        //        IXX += node.Properties.Ixx;
        //        IYY += node.Properties.Iyy;
        //        IZZ += node.Properties.Izz;
        //        IXY += node.Properties.Ixy;
        //        IYZ += node.Properties.Iyz;
        //        IZX += node.Properties.Izx;
        //    });
        //    if (Math.Abs(M) <= double.Epsilon)
        //    {
        //        this.ResetValue();
        //    }
        //    else
        //    {
        //        double mass = M;
        //        double xcg = MX / M;
        //        double ycg = MY / M;
        //        double zcg = MZ / M;
        //        double xxcg = xcg * xcg;
        //        double yycg = ycg * ycg;
        //        double zzcg = zcg * zcg;
        //        double xycg = xcg * ycg;
        //        double yzcg = ycg * zcg;
        //        double zxcg = zcg * xcg;
        //        double ixx = IXX + MYY + MZZ + (yycg + zzcg) * M - 2 * ycg * MY - 2 * zcg * MZ;
        //        double iyy = IYY + MZZ + MXX + (zzcg + xxcg) * M - 2 * zcg * MZ - 2 * xcg * MX;
        //        double izz = IZZ + MXX + MYY + (xxcg + yycg) * M - 2 * xcg * MX - 2 * ycg * MY;
        //        double ixy = IXY + xycg * M - xcg * MY - ycg * MX + MXY;
        //        double iyz = IYZ + yzcg * M - ycg * MZ - zcg * MY + MYZ;
        //        double izx = IZX + zxcg * M - zcg * MX - xcg * MZ + MZX;
        //        this.SetValue(mass, xcg, ycg, zcg, ixx, iyy, izz, ixy, iyz, izx);
        //    }
        //    this._updated = true;
        //}
    }

    public class MassProperties
    {
        public double Mass = 0;
        public double Xcg = 0;
        public double Ycg = 0;
        public double Zcg = 0;
        public double Ixx = 0;
        public double Iyy = 0;
        public double Izz = 0;
        public double Ixy = 0;
        public double Iyz = 0;
        public double Izx = 0;

        public void ResetValue()
        {
            this.Mass = 0;
            this.Xcg = 0;
            this.Ycg = 0;
            this.Zcg = 0;
            this.Ixx = 0;
            this.Iyy = 0;
            this.Ixy = 0;
            this.Iyz = 0;
            this.Izx = 0;
        }

        public void SetValue(double mass = double.NaN, 
            double xcg = double.NaN, double ycg = double.NaN, double zcg = double.NaN,
            double ixx = double.NaN, double iyy = double.NaN, double izz = double.NaN, 
            double ixy = double.NaN, double iyz = double.NaN, double izx = double.NaN)
        {
            if (mass != double.NaN)
                this.Mass = mass;
            if (xcg != double.NaN)
                this.Xcg = xcg;
            if (ycg != double.NaN)
                this.Ycg = ycg;
            if (zcg != double.NaN)
                this.Zcg = zcg;
            if (ixx != double.NaN)
                this.Ixx = ixx;
            if (iyy != double.NaN)
                this.Iyy = iyy;
            if (izz != double.NaN)
                this.Izz = izz;
            if (ixy != double.NaN)
                this.Ixy = ixy;
            if (iyz != double.NaN)
                this.Iyz = iyz;
            if (izx != double.NaN)
                this.Izx = izx;
        }

        public MassProperties(double mass,
            double xcg, double ycg, double zcg,
            double ixx, double iyy, double izz,
            double ixy, double iyz, double izx)
        {
            SetValue(mass, xcg, ycg, zcg, ixx, iyy, izz, ixy, iyz, izx);
        }
    }

    public class MassNodeForSample : MassNode
    {
        //用于对重量随机采样
        private RandModel _sampleModel;
        public RandModel SampleModel
        {
            get { return _sampleModel; }
        }
        public MassNodeForSample(double p = 1, double mu = 0, double sigma = 1)
        {
            this._sampleModel = new MixedModel(p, mu, sigma);
        }

        
        public List<MassProperties> SampleList = new List<MassProperties>();

        //public void Sample(int n)
        //{
        //    SampleList.Clear();
        //    if (this.ChildNodes.IsEmpty)
        //    {  
        //        Parallel.For(0, n, (index) =>
        //        {
        //            SampleList.Add(new MassProperties(SampleModel.Sample(), Xcg, Ycg, Zcg, Ixx, Iyy, Izz, Ixy, Iyz, Izx));
        //        });
        //        return;
        //    }
        //    Parallel.ForEach(ChildNodes, (node) =>
        //    {
        //        MassNodeForSample snode = node as MassNodeForSample;
        //        snode.Sample(n);
        //    });
        //    SampleList = new List<MassProperties>(n);
        //    Parallel.ForEach(SampleList, (sample) =>
        //    {
        //        sample = new MassProperties(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        //        int i = SampleList.IndexOf(sample);
        //        double M = 0;
        //        double MX = 0, MY = 0, MZ = 0;
        //        double MXX = 0, MYY = 0, MZZ = 0, MXY = 0, MYZ = 0, MZX = 0;
        //        double IXX = 0, IYY = 0, IZZ = 0, IXY = 0, IYZ = 0, IZX = 0;
        //        Parallel.ForEach(ChildNodes, (node) =>
        //        {
        //            MassNodeForSample snode = node as MassNodeForSample;
        //            double Mi = snode.SampleList[i].Mass;
        //            double MXi = Mi * snode.SampleList[i].Xcg;
        //            double MYi = Mi * snode.SampleList[i].Ycg;
        //            double MZi = Mi * snode.SampleList[i].Zcg;
        //            double MXXi = MXi * snode.SampleList[i].Xcg;
        //            double MYYi = MYi * snode.SampleList[i].Ycg;
        //            double MZZi = MZi * snode.SampleList[i].Zcg;
        //            double MXYi = MXi * snode.SampleList[i].Ycg;
        //            double MYZi = MYi * snode.SampleList[i].Zcg;
        //            double MZXi = MZi * snode.SampleList[i].Xcg;
        //            M += Mi;
        //            MX += MXi;
        //            MY += MYi;
        //            MZ += MZi;
        //            MXX += MXXi;
        //            MYY += MYYi;
        //            MZZ += MZZi;
        //            MXY += MXYi;
        //            MYZ += MYZi;
        //            MZX += MZXi;
        //            IXX += snode.SampleList[i].Ixx;
        //            IYY += snode.SampleList[i].Iyy;
        //            IZZ += snode.SampleList[i].Izz;
        //            IXY += snode.SampleList[i].Ixy;
        //            IYZ += snode.SampleList[i].Iyz;
        //            IZX += snode.SampleList[i].Izx;
        //        });
        //        if (Math.Abs(M) <= double.Epsilon)
        //        {
        //            sample.ResetValue();
        //        }
        //        else
        //        {
        //            double mass = M;
        //            double xcg = MX / M;
        //            double ycg = MY / M;
        //            double zcg = MZ / M;
        //            double xxcg = xcg * xcg;
        //            double yycg = ycg * ycg;
        //            double zzcg = zcg * zcg;
        //            double xycg = xcg * ycg;
        //            double yzcg = ycg * zcg;
        //            double zxcg = zcg * xcg;
        //            double ixx = IXX + MYY + MZZ + (yycg + zzcg) * M - 2 * ycg * MY - 2 * zcg * MZ;
        //            double iyy = IYY + MZZ + MXX + (zzcg + xxcg) * M - 2 * zcg * MZ - 2 * xcg * MX;
        //            double izz = IZZ + MXX + MYY + (xxcg + yycg) * M - 2 * xcg * MX - 2 * ycg * MY;
        //            double ixy = IXY + xycg * M - xcg * MY - ycg * MX + MXY;
        //            double iyz = IYZ + yzcg * M - ycg * MZ - zcg * MY + MYZ;
        //            double izx = IZX + zxcg * M - zcg * MX - xcg * MZ + MZX;
        //            sample.SetValue(mass, xcg, ycg, zcg, ixx, iyy, izz, ixy, iyz, izx);
        //        }
        //    });

        //}
    }
}
