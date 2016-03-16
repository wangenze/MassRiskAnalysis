using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MassRiskAnalysis
{
    public static class Numerics
    {
        static Random random = new Random();
        const double TWO_PI = 2.0 * Math.PI;
        public static double GenerateNormalSample(double mu, double sigma)
        {
            double z0;
            //double z1;
            double u1, u2;
            do
            {
                u1 = random.NextDouble();
                u2 = random.NextDouble();
            }
            while (u1 <= double.Epsilon);

            z0 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(TWO_PI * u2);
            //z1 = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(two_pi * u2);
            return z0 * sigma + mu;
        }

        public static double GenerateUniformSample(double a = 0, double b = 1)
        {
            double length = b - a;
            double u1 = random.NextDouble();
            return length * u1 - a;
        }
    }

    public abstract class RandModel
    {
        public abstract double Sample();
    }

    public class NormalModel : RandModel
    {
        private double _mu;
        public double Mu
        {
            get { return _mu; }
            set
            {
                if (!double.IsNaN(value))
                    this._mu = value;
                else
                    this._mu = 0;
            }
        }
        private double _sigma;
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (!double.IsNaN(value))
                    this._sigma = value;
                else
                    this._sigma = 0;
            }
        }
        public NormalModel(double mu = 0, double sigma = 1)
        {
            if (!double.IsNaN(mu))
                this._mu = mu;
            if (!double.IsNaN(sigma))
                this._sigma = sigma;
            if (this._sigma < 0.0) this._sigma = -this._sigma;
        }

        public override double Sample()
        {
            return Numerics.GenerateNormalSample(this._mu, this._sigma);
        }
    }

    public class UniformModel : RandModel
    {
        private double _a;
        public double A
        {
            get { return _a; }
            set
            {
                if (!double.IsNaN(value))
                    this._a = value;
                else
                    this._a = 0;
            }
        }
        private double _b;
        public double B
        {
            get { return _b; }
            set
            {
                if (!double.IsNaN(value))
                    this._b = value;
                else
                    this._b = 1;
            }
        }
        public UniformModel(double a = 0, double b = 1)
        {
            if (!double.IsNaN(a))
                this._a = a;
            if (!double.IsNaN(b))
                this._b = b;
        }

        public override double Sample()
        {
            return Numerics.GenerateUniformSample(_a, _b);
        }
    }

    public class MixedModel : RandModel
    {
        private double _p;
        public double P
        {
            get { return _p; }
            set
            {
                if ((!double.IsNaN(value)) && (value >= 0.0 && value <= 1.0))
                    this._p = value;
                else
                    this._p = 1;
            }
        }
        private double _mu;
        public double Mu
        {
            get { return _mu; }
            set
            {
                if (!double.IsNaN(value))
                    this._mu = value;
                else
                    this._mu = 0;
            }
        }
        private double _sigma;
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (!double.IsNaN(value))
                    this._sigma = value;
                else
                    this._sigma = 0;
            }
        }
        public MixedModel(double p = 1, double mu = 0, double sigma = 1)
        {
            if ((!double.IsNaN(p)) && (p >= 0.0 && p <= 1.0))
                this._p = p;
            if (!double.IsNaN(mu))
                this._mu = mu;
            if (!double.IsNaN(sigma))
                this._sigma = sigma;
            if (this._sigma < 0.0) this._sigma = -this._sigma;
        }
        public override double Sample()
        {
            if (this._p >= 1 - double.Epsilon || this._p >= Numerics.GenerateUniformSample())
                return Numerics.GenerateNormalSample(this._mu, this._sigma);
            else
                return 0.0;
        }
    }
}
