use fibonacci::nth_fibonacci;

fn main() {
    let n = 1_000_000;
    let nth_fib = nth_fibonacci(n);
    println!("The {}th Fibonacci number is: {}", n, nth_fib);
}
