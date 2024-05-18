/// Functions for efficiently working with Fibonacci numbers
/// 
/// The emphasis of the functions are to work on large
/// (numbers that won't fit into u64 types) Fibonacci numbers,
/// due to this constraint a BigUint implementation is used.
/// This will lead to inefficiency for numbers that could fit
/// into an u64 type.


use num_bigint::BigUint;
use num_traits::{One, Zero};

/// Multiply lhs by rhs
/// 
/// Both lhs and rhs must be flattened, row major, 2x2 matrices
fn matrix_multiply_2x2(lhs: &[BigUint; 4], rhs: &[BigUint; 4]) -> [BigUint; 4] {
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
            // Here we use exponentiation by squaring for 2x2 matrices
            // https://en.wikipedia.org/wiki/Exponentiation_by_squaring
            // TODO: How does Rust handle recursion?
            // TODO: Would this be faster as an iterative algorithm?
            let matrix_squared = matrix_multiply_2x2(&matrix, &matrix);
            let n_div_2 = n / 2;
            if n % 2 == 0 {
                matrix_power_2x2(&matrix_squared, n_div_2)
            } else {
                matrix_multiply_2x2(&matrix, &matrix_power_2x2(&matrix_squared, n_div_2))
            }
        }
    }
}

/// # Compute the nth Fibonacci number
///
/// This function uses the matrix power approach to calculate the nth Fibonacci number
/// in O(log n) time complexity. The function allows for the computation of very large
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
            let fib_matrix = matrix_power_2x2(&fibonacci_base, n - 1);
            fib_matrix[0].clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigUint;

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
}
