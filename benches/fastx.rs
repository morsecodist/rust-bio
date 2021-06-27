use bio::io::{fasta, fastx};
use rand::distributions::Alphanumeric;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::io;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

const BASES: &[u8] = b"ACTG";
const ID_LEN: usize = 10;
const DESC_LEN: usize = 20;
const SEQ_LEN: usize = 100;
const FASTA_SIZE: usize = 1000;

fn gen_random_fasta() -> io::Result<Vec<u8>> {
    let mut raw_writer = Vec::new();
    {
        let mut w = fasta::Writer::new(&mut raw_writer);
        for _ in 0..FASTA_SIZE {
            let id: String = StdRng::seed_from_u64(42)
                .sample_iter(&Alphanumeric)
                .take(ID_LEN)
                .map(char::from)
                .collect();

            let desc: String = StdRng::seed_from_u64(4)
                .sample_iter(&Alphanumeric)
                .take(DESC_LEN)
                .map(char::from)
                .collect();

            let seq: Vec<u8> = (0..SEQ_LEN)
                .map(|_| {
                    let idx = StdRng::seed_from_u64(65).gen_range(0..BASES.len());
                    BASES[idx]
                })
                .collect();

            w.write_record(&fasta::Record::with_attrs(&id, Some(&desc), &seq))?;
        }
    }
    Ok(raw_writer)
}

fn fastx_count_bases<T, E, I>(records: I) -> Result<usize, E>
where
    T: fastx::Record,
    E: std::error::Error,
    I: fastx::Records<T, E>,
{
    let mut nb_bases = 0;
    for result in records {
        let record = result?;
        nb_bases += record.seq().len();
    }
    Ok(nb_bases)
}

fn fasta_count_bases<R>(records: fasta::Records<R>) -> io::Result<usize>
where
    R: io::Read,
{
    let mut nb_bases = 0;
    for result in records {
        let record = result?;
        nb_bases += record.seq().len();
    }
    Ok(nb_bases)
}

fn fastx_check<T, E, I>(records: I) -> Result<(), String>
where
    T: fastx::Record,
    E: std::error::Error,
    I: fastx::Records<T, E>,
{
    for result in records {
        let record = result.map_err(|e| format!("{}", e))?;
        record.check().map_err(|e| e.to_owned())?;
    }
    Ok(())
}

fn fasta_check<R>(records: fasta::Records<R>) -> Result<(), String>
where
    R: io::Read,
{
    for result in records {
        let record = result.map_err(|e| format!("{}", e))?;
        record.check().map_err(|e| e.to_owned())?;
    }
    Ok(())
}

fn bench_fasta_count(c: &mut Criterion) {
    let mut data = io::Cursor::new(gen_random_fasta().unwrap());
    c.bench_function("fasta count bases", |b| {
        b.iter(|| {
            data.set_position(0);
            let records = fasta::Reader::new(&mut data).records();
            fasta_count_bases(records).unwrap()
        })
    });
}

fn bench_fastx_fasta_count(c: &mut Criterion) {
    let mut data = io::Cursor::new(gen_random_fasta().unwrap());
    c.bench_function("fastx fasta count", |b| {
        b.iter(|| {
            data.set_position(0);
            let records = fasta::Reader::new(&mut data).records();
            fastx_count_bases(records).unwrap();
        })
    });
}

fn bench_either_fasta_count(c: &mut Criterion) {
    let mut data = io::Cursor::new(gen_random_fasta().unwrap());
    c.bench_function("either fasta count", |b| {
        b.iter(|| {
            data.set_position(0);
            let records = fastx::EitherRecords::new(&mut data);
            fastx_count_bases(records).unwrap();
        })
    });
}

fn bench_fasta_check(c: &mut Criterion) {
    let mut data = io::Cursor::new(gen_random_fasta().unwrap());
    c.bench_function("fasta check", |b| {
        b.iter(|| {
            data.set_position(0);
            let records = fasta::Reader::new(&mut data).records();
            fasta_check(records).unwrap();
        })
    });
}

fn bench_fastx_fasta_check(c: &mut Criterion) {
    let mut data = io::Cursor::new(gen_random_fasta().unwrap());
    c.bench_function("fastx fasta check", |b| {
        b.iter(|| {
            data.set_position(0);
            let records = fasta::Reader::new(&mut data).records();
            fastx_check(records).unwrap();
        })
    });
}

fn bench_either_fasta_check(c: &mut Criterion) {
    let mut data = io::Cursor::new(gen_random_fasta().unwrap());
    c.bench_function("either fasta check", |b| {
        b.iter(|| {
            data.set_position(0);
            let records = fastx::EitherRecords::new(&mut data);
            fastx_check(records).unwrap();
        })
    });
}

criterion_group!(benches, bench_fasta_count, bench_fastx_fasta_count, bench_either_fasta_count, bench_fasta_check, bench_fastx_fasta_check, bench_either_fasta_check);
criterion_main!(benches);
