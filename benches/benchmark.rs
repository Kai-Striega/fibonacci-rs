use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fibonacci::nth_fibonacci;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("fib 1_000_000", |b| {
        b.iter(|| nth_fibonacci(black_box(1_000_000)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
