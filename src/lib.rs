/// Functions for efficiently working with Fibonacci numbers
///
/// The Fibonacci series is the series of numbers where each number
/// is the sum of the two preceding numbers, starting with:
///
/// F_0, F_1 = 0, 1
///
/// This gives the sequence 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, ...
///
/// ## Warning:
///
/// The emphasis of the functions are to work on large
/// (numbers that won't fit into u64 type) Fibonacci numbers,
/// due to this constraint a BigUint implementation is used.
/// This will lead to inefficiency for numbers that could fit
/// into an u64 type.
use num_bigint::{BigUint, ToBigUint};
use num_traits::{One, Zero};

/// Multiply lhs by rhs
///
/// Both lhs and rhs must be flattened, row major, 2x2 matrices
fn matrix_multiply_2x2(lhs: &[BigUint; 4], rhs: &[BigUint; 4]) -> [BigUint; 4] {
    // TODO: Optimise?
    //       For Fibonacci matrices ``a12`` and ``a21`` will always be the same,
    //       we can calculate this once and use the same result.
    let a11 = &lhs[0] * &rhs[0] + &lhs[1] * &rhs[2];
    let a12 = &lhs[0] * &rhs[1] + &lhs[1] * &rhs[3];
    let a21 = &lhs[2] * &rhs[0] + &lhs[3] * &rhs[2];
    let a22 = &lhs[2] * &rhs[1] + &lhs[3] * &rhs[3];
    [a11, a12, a21, a22]
}

/// Raise ``matrix`` to the power of ``n``
///
/// The matrix must be a flattened, row major 2x2 matrix
fn matrix_power_2x2(matrix: &[BigUint; 4], n: usize) -> [BigUint; 4] {
    match n {
        // any matrix to the power of 0 is the identity matrix.
        0 => [
            BigUint::one(),
            BigUint::zero(),
            BigUint::zero(),
            BigUint::one(),
        ],
        // any matrix to the power of 1 is itself.
        1 => matrix.clone(),
        _ => {
            // We use exponentiation by squaring for 2x2 matrices
            // https://en.wikipedia.org/wiki/Exponentiation_by_squaring
            let mut z = matrix.clone();
            let mut result = matrix.clone();
            let mut i = n - 1;

            while i > 0 {
                // i divmod 2, can this be done in one op?
                let bit = i % 2;
                i /= 2;
                
                if bit == 1 {
                    result = matrix_multiply_2x2(&result, &z);
                }

                z = matrix_multiply_2x2(&z, &z);
            }
            result
        }
    }
}

/// # Compute the nth Fibonacci number
///
/// ``nth_fibonacci`` uses the matrix power approach to calculate the nth Fibonacci number
/// in O(log n) time. ``nth_fibonacci`` allows for the computation of arbitrarily large
/// Fibonacci numbers due to using arbitrary sized unsigned integers. This comes with
/// a slight performance penalty for Fibonacci numbers that could fit into an u64 type.
///
/// # Example
///
/// ```
/// use num_bigint::ToBigUint;
/// use fibonacci::nth_fibonacci;
/// let n = 80;
/// let nth_fib = nth_fibonacci(n);
/// assert_eq!(nth_fib, 23_416_728_348_467_685u64.to_biguint().unwrap())
/// ```
pub fn nth_fibonacci(n: usize) -> BigUint {
    match n {
        0 => BigUint::zero(),
        _ => {
            let fibonacci_base = [
                BigUint::one(),
                BigUint::one(),
                BigUint::one(),
                BigUint::zero(),
            ];
            // We only need an exponent of n - 1 as the matrix contains F_{n+1}
            let exp = n - 1;
            let fib_matrix = matrix_power_2x2(&fibonacci_base, exp);
            fib_matrix[0].clone()
        }
    }
}

/// Check if a number is a perfect square
///
/// A perfect square is an element that is equal to the square of another element
fn is_perfect_square(n: &BigUint) -> bool {
    // TODO: Why is * needed before the rhs but not lhs?
    n.sqrt() * n.sqrt() == *n
}

/// Check if a number is a Fibonacci number
///
/// # Examples
/// ```
/// use num_bigint::ToBigUint;
/// use fibonacci::is_fibonacci;
///
/// let n = 55.to_biguint().unwrap();
/// let is_fib = is_fibonacci(&n);
/// assert_eq!(is_fib, true);
///
/// let n = 56.to_biguint().unwrap();
/// let is_fib = is_fibonacci(&n);
/// assert_eq!(is_fib, false);
/// ```
pub fn is_fibonacci(n: &BigUint) -> bool {
    // TODO should these be static?
    let zero = BigUint::ZERO;
    let four = 4.to_biguint().unwrap();
    let five = 5.to_biguint().unwrap();

    if n == &zero {
        return true;
    }

    // I don't like defining lots of variables like this,
    // but couldn't think of a neat way of doing it.
    let five_n_squared = five * n * n;
    let five_n_square_p_four = &five_n_squared + &four;
    let five_n_square_m_four = &five_n_squared - &four;
    is_perfect_square(&five_n_square_p_four) || is_perfect_square(&five_n_square_m_four)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_fib_0() {
        let result = nth_fibonacci(0);
        assert_eq!(result, BigUint::zero());
    }

    #[test]
    fn test_small_fib_1() {
        let result = nth_fibonacci(1);
        assert_eq!(result, BigUint::one());
    }

    #[test]
    fn test_small_fib_2() {
        let result = nth_fibonacci(2);
        assert_eq!(result, BigUint::one());
    }

    #[test]
    fn test_small_fib_3() {
        let result = nth_fibonacci(3);
        assert_eq!(result, 2.to_biguint().unwrap());
    }

    #[test]
    fn test_large_fib_80() {
        let result = nth_fibonacci(80);
        assert_eq!(result, 23_416_728_348_467_685u64.to_biguint().unwrap());
    }

    #[test]
    fn test_zero_is_fib() {
        let n = BigUint::zero();
        let is_fib = is_fibonacci(&n);
        assert_eq!(is_fib, true)
    }

    #[test]
    fn test_one_is_fib() {
        let n = BigUint::one();
        let is_fib = is_fibonacci(&n);
        assert_eq!(is_fib, true)
    }

    #[test]
    fn test_large_number_is_fib() {
        let n = nth_fibonacci(80);
        let is_fib = is_fibonacci(&n);
        assert_eq!(is_fib, true)
    }

    #[test]
    fn test_large_number_is_not_fib() {
        let n = nth_fibonacci(80) + BigUint::one();
        let is_fib = is_fibonacci(&n);
        assert_eq!(is_fib, false)
    }
}
