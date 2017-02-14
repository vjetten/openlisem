#ifndef GL3DMATH_H
#define GL3DMATH_H

#include <QVector3D>
#include <QVector2D>

class Vertex
{
public:
  // Constructors
  Q_DECL_CONSTEXPR Vertex();
  Q_DECL_CONSTEXPR explicit Vertex(const QVector3D &position);
  Q_DECL_CONSTEXPR Vertex(const QVector3D &position, const QVector3D &color);
  Q_DECL_CONSTEXPR Vertex(const QVector3D &position, const QVector3D &color, const QVector3D &normal);

  // Accessors / Mutators
  Q_DECL_CONSTEXPR const QVector3D& position() const;
  Q_DECL_CONSTEXPR const QVector3D& color() const;
  Q_DECL_CONSTEXPR const QVector3D& normal() const;
  void setPosition(const QVector3D& position);
  void setColor(const QVector3D& color);
  void setNormal(const QVector3D& normal);

  // OpenGL Helpers
  static const int PositionTupleSize = 3;
  static const int ColorTupleSize = 3;
  static const int NormalTupleSize = 3;
  static Q_DECL_CONSTEXPR int positionOffset();
  static Q_DECL_CONSTEXPR int colorOffset();
  static Q_DECL_CONSTEXPR int normalOffset();
  static Q_DECL_CONSTEXPR int stride();


  QVector3D m_position;
  QVector2D m_UV;
  QVector3D m_normal;

  QVector3D m_tangent;
  QVector3D m_bitangent;

  QVector3D m_color;

private:


};

/*******************************************************************************
 * Inline Implementation
 ******************************************************************************/

// Note: Q_MOVABLE_TYPE means it can be memcpy'd.
Q_DECLARE_TYPEINFO(Vertex, Q_MOVABLE_TYPE);

// Constructors
Q_DECL_CONSTEXPR inline Vertex::Vertex() {}
Q_DECL_CONSTEXPR inline Vertex::Vertex(const QVector3D &position) : m_position(position) {}
Q_DECL_CONSTEXPR inline Vertex::Vertex(const QVector3D &position, const QVector3D &color) : m_position(position), m_color(color) {}
Q_DECL_CONSTEXPR inline Vertex::Vertex(const QVector3D &position, const QVector3D &color, const QVector3D &normal) : m_position(position), m_color(color), m_normal(normal) {}

// Accessors / Mutators
Q_DECL_CONSTEXPR inline const QVector3D& Vertex::position() const { return m_position; }
Q_DECL_CONSTEXPR inline const QVector3D& Vertex::color() const { return m_color; }
Q_DECL_CONSTEXPR inline const QVector3D& Vertex::normal() const { return m_normal; }
void inline Vertex::setPosition(const QVector3D& position) { m_position = position; }
void inline Vertex::setColor(const QVector3D& color) { m_color = color; }
void inline Vertex::setNormal(const QVector3D& normal) { m_normal = normal; }

// OpenGL Helpers
Q_DECL_CONSTEXPR inline int Vertex::positionOffset() { return offsetof(Vertex, m_position); }
Q_DECL_CONSTEXPR inline int Vertex::colorOffset() { return offsetof(Vertex, m_color); }
Q_DECL_CONSTEXPR inline int Vertex::normalOffset() { return offsetof(Vertex, m_normal); }
Q_DECL_CONSTEXPR inline int Vertex::stride() { return sizeof(Vertex); }

class Vertex2D
{
public:
  // Constructors
  Q_DECL_CONSTEXPR Vertex2D();
  Q_DECL_CONSTEXPR explicit Vertex2D(const QVector2D &position);

  // Accessors / Mutators
  Q_DECL_CONSTEXPR const QVector2D& position() const;

  void setPosition(const QVector2D& position);

  // OpenGL Helpers
  static const int PositionTupleSize = 2;
  static Q_DECL_CONSTEXPR int positionOffset();
  static Q_DECL_CONSTEXPR int stride();

private:
  QVector2D m_position;

};

/*******************************************************************************
 * Inline Implementation
 ******************************************************************************/

// Note: Q_MOVABLE_TYPE means it can be memcpy'd.
Q_DECLARE_TYPEINFO(Vertex2D, Q_MOVABLE_TYPE);

// Constructors
Q_DECL_CONSTEXPR inline Vertex2D::Vertex2D() {}
Q_DECL_CONSTEXPR inline Vertex2D::Vertex2D(const QVector2D &position) : m_position(position) {}

// Accessors / Mutators
Q_DECL_CONSTEXPR inline const QVector2D& Vertex2D::position() const { return m_position; }
void inline Vertex2D::setPosition(const QVector2D& position) { m_position = position; }
// OpenGL Helpers
Q_DECL_CONSTEXPR inline int Vertex2D::positionOffset() { return offsetof(Vertex2D, m_position); }
Q_DECL_CONSTEXPR inline int Vertex2D::stride() { return sizeof(Vertex2D); }


// THIS CLASS IS A TRANSLATION TO C++11 FROM THE REFERENCE
// JAVA IMPLEMENTATION OF THE IMPROVED PERLIN FUNCTION (see http://mrl.nyu.edu/~perlin/noise/)
// THE ORIGINAL JAVA IMPLEMENTATION IS COPYRIGHT 2002 KEN PERLIN

// I (https://www.solarianprogrammer.com/2012/07/18/perlin-noise-cpp-11/) ADDED AN EXTRA METHOD THAT GENERATES A NEW PERMUTATION VECTOR (THIS IS NOT PRESENT IN THE ORIGINAL IMPLEMENTATION)

class PerlinNoise {
    // The permutation vector
    std::vector<int> p;
public:
    // Initialize with the reference values for the permutation vector
    PerlinNoise();
    // Generate a new permutation vector based on the value of seed
    PerlinNoise(unsigned int seed);
    // Get a noise value, for 2D images z can have any value
    double noise(double x, double y, double z);
private:
    double fade(double t);
    double lerp(double t, double a, double b);
    double grad(int hash, double x, double y, double z);
};


#endif // GL3DMATH_H


