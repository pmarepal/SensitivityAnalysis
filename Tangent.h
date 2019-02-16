#ifndef _TANGENT_H_
#define _TANGENT_H_


class Tangent 
{
public:

  typedef Tangent This_T;
  typedef Tangent T_Scalar;

  
  // no initialization for empty ctor - consistent with builtins
  Tangent():
  _v(0),
  _dv1(0),
  _dv2(0),
  _dv1dv2(0)
  {}
  
  
  Tangent(const double v, const double dv1, const double dv2, const double dv1dv2) :
    _v(v),
    _dv1(dv1),
    _dv2(dv2),
    _dv1dv2(dv1dv2)
  {}

  Tangent(const double v) :
    _v(v),
    _dv1(0.),
    _dv2(0.),
    _dv1dv2(0)
  {}


  //to allow statements like T x(0)
  Tangent(const int v) :
    _v(double(v)),
    _dv1(0.),
    _dv2(0.),
    _dv1dv2(0.)
  {}

  Tangent(const Tangent& o) :
    _v(o._v),
    _dv1(o._dv1),
    _dv2(o._dv2),
    _dv1dv2(o._dv1dv2)
  {}
  
  ~Tangent()
  {}

  
  Tangent& operator=(const Tangent& o)
  {
    if (this == &o)
      return *this;
    _v = o._v;
    _dv1 = o._dv1;
    _dv2 = o._dv2;
    _dv1dv2 = o._dv1dv2;
    return *this;
  }

  Tangent& operator=(const double& f)
  {
    _v = f;
    _dv1 = 0.;
    _dv2 = 0.;
    _dv1dv2 = 0.;
    return *this;
  }

  // to allow statements like x = 0;
  Tangent& operator=(const int& i)
  {
    _v = double(i);
    _dv1 = 0;
    _dv2 = 0;
    _dv1dv2 = 0;
    return *this;
  }

  Tangent& operator+=(const Tangent& o)
  {
    _v += o._v;
    _dv1 += o._dv1;
    _dv2 += o._dv2;
    _dv1dv2 += o._dv1dv2;
    return *this;
  }

  Tangent& operator-=(const Tangent& o)
  {
    _v -= o._v;
    _dv1 -= o._dv1;
    _dv2 -= o._dv2;
    _dv1dv2 -= o._dv1dv2;
    return *this;
  }

  Tangent& operator*=(const Tangent& o)
  {
    _dv1dv2 = _dv1dv2*o._v + o._dv1dv2*_v + _dv1*o._dv2 + o._dv1*_dv2;
    _dv1 = _dv1*o._v + _v*o._dv1;
    _dv2 = _dv2*o._v + _v*o._dv2;
    _v *= o._v;
    return *this;
  }

  Tangent& operator/=(const Tangent& o)
  {
    _dv1dv2 = (_dv1dv2*o._v*o._v*o._v - o._dv1dv2*_v*o._v*o._v - _dv1*o._dv2*o._v*o._v - o._dv1*_dv2*o._v*o._v + 2*o._dv1*o._dv2*_v*o._v)/(o._v*o._v*o._v*o._v);
    _dv1 = (o._v*_dv1 - _v*o._dv1)/(o._v*o._v);
    _dv2 = (o._v*_dv2 - _v*o._dv2)/(o._v*o._v);
    _v /= o._v;
    return *this;
  }


  Tangent& operator+=(const double& o)
  {
    _v += o;
    return *this;
  }

  Tangent& operator-=(const double& o)
  {
    _v -= o;
    return *this;
  }

  Tangent& operator*=(const double& o)
  {
    _v *= o;
    _dv1 *= o;
    _dv2 *= o; //check this
    _dv1dv2 *= o;
    return *this;
  }

  Tangent& operator/=(const double& o)
  {
    _v /= o;
    _dv1 /= o;
    _dv2 /= o;
    _dv1dv2 /= o;
    return *this;
  }

  Tangent& operator+=(const int& i)
  {
    _v += double(i);
    return *this;
  }
  
  Tangent& operator-=(const int& i)
  {
    _v -= double(i);
    return *this;
  }
  
  Tangent& operator*=(const int& i)
  {
    _v *= double(i);
    _dv1 *= double(i); //check this too
    _dv2 *= double(i);
    _dv1dv2 *= double(i); 
    return *this;
  }
  
  Tangent& operator/=(const int& i)
  {
    _v /= double(i);
    _dv1 /= double(i); //check this
    _dv2 /= double(i);
    _dv1dv2 /= double(i);
    return *this;
  }

#define TANGENT_RELATIONAL_OPS(opname,_op_)     \
  bool opname(const Tangent& o) const           \
  {                                             \
    return (_v _op_ o._v);                      \
  }                                             \
    bool opname(const double& o) const          \
  {                                             \
    return (_v _op_ o);                         \
  }                                             \
    bool opname(const int& o) const             \
  {                                             \
    return (_v _op_ double(o));                 \
  }

  TANGENT_RELATIONAL_OPS(operator>,>);
  TANGENT_RELATIONAL_OPS(operator>=,>=);
  TANGENT_RELATIONAL_OPS(operator<,<);
  TANGENT_RELATIONAL_OPS(operator<=,<=);
  TANGENT_RELATIONAL_OPS(operator==,==);
  TANGENT_RELATIONAL_OPS(operator!=,!=);
  
#undef TANGENT_RELATIONAL_OPS
  
  void print(std::ostream &os) const
  {
    os << "< "  << _v << " , " << _dv1 <<  " , " << _dv2 << " , " << _dv1dv2 << "> ";
  }

  static Tangent getZero()
  {
    double zero = 0;
    return Tangent(zero,zero,zero,zero);
  }

  static Tangent getUnity()
  {
    return Tangent(1,0,0,0);
  }

  static Tangent getNegativeUnity()
  {
    return Tangent(-1,0,0,0);
  }


  Tangent fabs() const
  {
    return Tangent(::fabs(_v), (_v > 0 ? _dv1 : -_dv1), (_v > 0 ? _dv2 : -_dv2), (_v > 0 ? _dv1dv2 : -_dv1dv2));
  }
   
  
  double _v;
  double _dv1;
  double _dv2;
  double _dv1dv2;
};


#define TANGENT_BINARY_OP(opname,_op_)                  \
  Tangent opname(const Tangent& a, const Tangent& b)    \
  {                                                     \
    return Tangent(a) _op_ b;                            \
  }                                                      \
  Tangent opname(const Tangent& a, const double& b)      \
  {                                                      \
    return Tangent(a) _op_ b;                            \
  }                                                      \
  Tangent opname(const double& a, const Tangent& b)      \
  {                                                      \
    return Tangent(a) _op_ b;                            \
  }                                                      \
  Tangent opname(const Tangent& a, const int& b)         \
  {                                                      \
    return Tangent(a) _op_ b;                            \
  }                                                      \
  Tangent opname(const int& a, const Tangent& b)         \
  {                                                      \
    return Tangent(a) _op_ b;                            \
  }                                                      \
                                                         
TANGENT_BINARY_OP(operator+,+=);
TANGENT_BINARY_OP(operator-,-=);
TANGENT_BINARY_OP(operator*,*=);
TANGENT_BINARY_OP(operator/,/=);

Tangent operator+(const Tangent& a)
{
  return a;
}

Tangent operator-(const Tangent& a)
{
  return Tangent(-a._v,-a._dv1,-a._dv2,-a._dv1dv2);
}

inline std::ostream &operator<<(std::ostream &os, const Tangent &a)
{
  a.print(os);
  return os;
}

Tangent sin(const Tangent& a)
{
  return Tangent(sin(a._v), a._dv1*cos(a._v), a._dv2*cos(a._v), a._dv1dv2*cos(a._v)+a._dv1*a._dv2*-sin(a._v));
}

Tangent cos(const Tangent& a)
{
 return Tangent(cos(a._v), a._dv1*-sin(a._v), a._dv2*-sin(a._v), -a._dv1dv2*sin(a._v)+a._dv1*a._dv2*-cos(a._v));
}

Tangent exp(const Tangent& a)
{
 return Tangent(exp(a._v), a._dv1*exp(a._v), a._dv2*exp(a._v), a._dv1dv2*exp(a._v)+a._dv1*a._dv2*exp(a._v));
}

Tangent pow(const Tangent& a, const double& b)
{
  double arg1 = a._v == 0.0 ? 0 : pow(a._v,b); 
  double arg2 = a._v == 0.0 ? 0 : a._dv1*b*pow(a._v,b-1);
  double arg3 = a._v == 0.0 ? 0 : a._dv2*b*pow(a._v,b-1);
  double arg4 = a._v == 0.0 ? 0 : a._dv1*a._dv2*b*(b-1)*pow(a._v,b-2)+a._dv1dv2*b*pow(a._v,b-1);
  return Tangent(arg1, arg2, arg3, arg4);
}

Tangent fabs(const Tangent& a)
{
   return Tangent(::fabs(a._v), (a._v > 0 ? a._dv1 : -a._dv1), (a._v > 0 ? a._dv2 : -a._dv2), (a._v > 0 ? a._dv1dv2 : -a._dv1dv2));
}
  
Tangent abs(const Tangent& a)
{
   return Tangent(::fabs(a._v), (a._v > 0 ? a._dv1 : -a._dv1), (a._v > 0 ? a._dv2 : -a._dv2), (a._v > 0 ? a._dv1dv2 : -a._dv1dv2));
}

Tangent sqrt(const Tangent& a)
{
//  double sqv = sqrt(a._v);
//  return Tangent(sqv, sqv==0.0 ? 0 : a._dv1*0.5/sqv, sqv==0.0 ? 0 : a._dv2*0.5/sqv, sqv==0.0 ? 0 : a._dv1dv2*0.5*pow(a._v,-0.5)+a._dv1*a._dv2*-1*0.25*pow(a._v,-1.5));
   Tangent sqv;
   sqv = pow(a,0.5);
   return sqv;
}


double ceil(const Tangent& a)
{
  return ceil(a._v);
}

Tangent log(const Tangent& a)
{
  // FIXME
  return Tangent(log(a._v),0,0,0);
}


inline const Tangent& conj(const Tangent& x) { return x; }
inline const Tangent& real(const Tangent& x) { return x; }
inline Tangent imag(const Tangent& x) { return Tangent(x._v,0.,0.,0.); }
//inline Tangent abs(const Tangent& x) { return fabs(x); }
inline Tangent abs2(const Tangent& x) { return x*x; }

namespace Eigen{

template<> struct NumTraits<Tangent>
: NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
typedef Tangent Real;
typedef Tangent NonInteger;
typedef Tangent Nested;
enum {
IsComplex = 0,
IsInteger = 0,
IsSigned = 1,
RequireInitialization = 1,
ReadCost = 1,
AddCost = 3,
MulCost = 3
};
};

namespace internal {

template<>
struct significant_decimals_default_impl<Tangent, false>
{
  typedef double RealScalar;
  static inline int run()
  {
    using std::ceil;
    return cast<RealScalar,int>(ceil(-log(NumTraits<RealScalar>::epsilon())/log(RealScalar(10))));
  }
};

template<>
struct significant_decimals_default_impl<std::complex<Tangent>, false>
{
  typedef double RealScalar;
  static inline int run()
  {
    using std::ceil;
    return cast<RealScalar,int>(ceil(-log(NumTraits<RealScalar>::epsilon())/log(RealScalar(10))));
  }
};

};
};

#endif
