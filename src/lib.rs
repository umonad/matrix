use std::ops::{Add, Mul, Sub}; //

// _____________________________________________________STRUCTURES_VECTOR_AND_MATRIX____________________________________________________________//

#[derive(Debug, Clone, PartialEq)]
pub struct Vector<const N: usize> {
    //TODO: replace f64 by field -> what is field
    pub data: [f64; N], // [type of vale, size of the array]
}
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<const W: usize, const H: usize> {
    //TODO: replace f64 by field -> what is field
    pub data: [[f64; W]; H], // [type of vale, size of the line of the matrice, size of the matrice]
}

// _____________________________________________________VECTOR IMPLEMENTATION____________________________________________________________//

impl<const N: usize> Vector<N> {
    // NEW VECTOR
    pub fn new(data: [f64; N]) -> Self {
        Self { data }
    }

    fn dot(self, u: Vector<N>) -> f64 {
        (self * u).data.iter().sum()
    }

    fn norm_1(self) -> f64 {
        self.data.map(f64::abs).iter().sum()
    }
    fn norm(self) -> f64 {
        (self.data.map(|x| x * x).iter().sum::<f64>()).sqrt()
    }
    fn norm_inf(self) -> f64{
        self.data.iter().copied().map(f64::abs).fold(0.0_f64, f64::max)
    }
}

//TODO: impl "&self macro" for stop clonning nd impl operator with both value and copie of value
impl<const N: usize> Add for Vector<N> {
    // VECTOR + VECTOR
    type Output = Vector<N>;

    fn add(self, other: Self) -> Self::Output {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = other.data[i] + self.data[i];
        }
        return Vector::new(res);
    }
}
//TODO: impl "&self macro" for stop clonning nd impl operator with both value and copie of value
impl<const N: usize> Sub for Vector<N> {
    // VECTOR - VECTOR
    type Output = Vector<N>;

    fn sub(self, other: Self) -> Self::Output {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = other.data[i] - self.data[i];
        }
        return Vector::new(res);
    }
}
impl<const N: usize> Mul<f64> for Vector<N> {
    // VECTOR * SCALAR
    type Output = Vector<N>;

    fn mul(self, other: f64) -> Self::Output {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = other * self.data[i];
        }
        return Vector::new(res);
    }
}

impl<const N: usize> Mul<Vector<N>> for f64 {
    type Output = Vector<N>;
    // SCALAR * VECTOR
    fn mul(self, other: Vector<N>) -> Self::Output {
        return other * self;
    }
}

//TODO: impl "&self macro" for stop clonning nd impl operator with both value and copie of value
impl<const N: usize> Mul for Vector<N> {
    // VECTOR * VECTOR
    type Output = Vector<N>;

    fn mul(self, other: Self) -> Self::Output {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = other.data[i] * self.data[i];
        }
        return Vector::new(res);
    }
}

fn linear_combination<const N: usize>(u: &[Vector<N>], coefs: &[f64]) -> Vector<N> {
    let mut res = [0.; N];
    for i in 0..u.len() {
        let mut sum = 0.;
        for j in 0..N {
            sum += coefs[i] * u[i].data[j];
        }
        res[i] = sum;
    }
    return Vector::new(res);
}

// _____________________________________________________MATRIX IMPLEMENTATION____________________________________________________________//

impl<const W: usize, const H: usize> Matrix<W, H> {
    // NEW MATRIX
    pub fn new(data: [[f64; W]; H]) -> Self {
        Self { data }
    }
}

// const  <> valeur connu a la compilation
impl<const W: usize, const H: usize> Add for Matrix<W, H> {
    // MATRIX + MATRIX
    type Output = Matrix<W, H>;

    fn add(self, matrice2: Matrix<W, H>) -> Self::Output {
        let mut res = [[0.0; W]; H];
        for i in 0..H {
            //lines
            for j in 0..W {
                //columns
                res[i][j] = self.data[i][j] + matrice2.data[i][j];
            }
        }
        return Matrix::new(res);
    }
}
// const  <> valeur connu a la compilation
impl<const W: usize, const H: usize> Sub for Matrix<W, H> {
    // MATRIX - MATRIX
    type Output = Matrix<W, H>;

    fn sub(self, matrice2: Matrix<W, H>) -> Self::Output {
        let mut res = [[0.0; W]; H];
        for i in 0..H {
            //lines
            for j in 0..W {
                //columns
                res[i][j] = self.data[i][j] - matrice2.data[i][j];
            }
        }
        return Matrix::new(res);
    }
}

impl<const W: usize, const H: usize> Mul<f64> for Matrix<W, H> {
    // MATRIX * SCALAR
    type Output = Matrix<W, H>;

    fn mul(self, scalar: f64) -> Self::Output {
        let mut res = [[0.0; W]; H];
        for i in 0..H {
            for j in 0..W {
                res[i][j] = self.data[i][j] * scalar
            }
        }
        return Matrix::new(res);
    }
}

impl<const W: usize, const H: usize> Mul<Matrix<W, H>> for f64 {
    type Output = Matrix<W, H>;
    // SCALAR * MATRIX
    fn mul(self, matrix: Matrix<W, H>) -> Self::Output {
        // delegate to the existing Matrix * f64 implementation
        matrix * self
    }
}

// const  <> valeur connu a la compilation
impl<const H1: usize, const W2: usize, const H2: usize> Mul<Matrix<W2, H2>> for Matrix<H2, H1> {
    // MATRIX * MATRIX
    type Output = Matrix<W2, H1>;

    fn mul(self, matrice2: Matrix<W2, H2>) -> Self::Output {
        let mut res = [[0.0; W2]; H1];
        for i in 0..H1 {
            //lines
            for j in 0..W2 {
                //columns
                let mut sum = 0.0;
                for k in 0..H2 {
                    // lines m1         //colomns m2
                    sum += self.data[i][k] * matrice2.data[k][j]; // la matrice d'application est self(m1), elle agit sur matrice2(m2)
                    // println!("{}", self.data[k][j]);
                }
                res[i][j] = sum
            }
        }
        return Matrix::new(res);
    }
}

fn lerp<V: Mul<f64, Output = V> + Add<V, Output = V>>(u: V, v: V, t: f64) -> V {
    u * (1. - t) + v * t
}

#[test]
fn testadd() {
    let v1 = Vector { data: [10.0, 20.] }; //without new method
    println!("{:?}", v1.clone() * 7.675);
    println!("{:?}", 1. * v1.clone());
    let v2 = Vector::new([10.0, 20.]); //with new method
    assert_eq!(v1.clone() + v2.clone(), Vector::new([20., 40.]));
    assert_eq!(v1.clone() - v2.clone(), Vector::new([0., 0.]));
    assert_eq!(v1.clone() * v2.clone(), Vector::new([100., 400.]));
    let m1 = Matrix::new([[0.0, -1.], [1.0, 0.]]);
    let m3 = Matrix::new([[5.0, -1.], [6.0, 8.]]);
    let m4 = Matrix::new([[7.0, -5.], [-9.0, 4.]]);
    let m2 = Matrix::new([[-3.0], [-2.]]);
    println!("{:?}", m1.clone() * m2.clone());
    assert_eq!(
        m3.clone() + m4.clone(),
        Matrix::new([[12., -6.], [-3., 12.]])
    );
    assert_eq!(m3.clone() * 1., m3.clone());
    assert_eq!(2. * m3.clone(), m3.clone() * 2.);
    println!("{:?}", m3.clone() * 2.);
}
