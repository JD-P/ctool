extern crate rand;

use structopt::StructOpt;
use std::collections::HashMap;
use std::collections::HashSet;
use regex::Regex;
use rand::{thread_rng, Rng};

const ALPHABET:[&str;26] = ["A","B","C",
                            "D","E","F",
                            "G","H","I",
                            "J","K","L",
                            "M","N","O",
                            "P","Q","R",
                            "S","T","U",
                            "V","W","X",
                            "Y","Z"];

/// Use frequency analysis to derive a key for ngram simple substitution ciphers
#[derive(StructOpt)]
#[structopt(name="Cipher Tool", about="Utility for analyzing cryptograms.")]
struct Cli {
    /// The ciphertext to analyze.
    cryptogram: String,
    /// The n-gram value for how many ciphertext characters correspond to one plaintext character.
    #[structopt(long="ngram-ssc", default_value="1")]
    ngram_ssc: usize,
    /// The n-gram value we use to analyze frequency for simple substitution, e.g trigram analysis is 3.
    #[structopt(short="n", long="ngram-freq", default_value="3")]
    n_gram: usize,
}

fn build_plaintext<'a>(ciphertext: &'a String, key: &[String;26]) -> String {
    let mut plaintext: String = "".to_string();
    let mut key_map = HashMap::new();
    let mut ngram: usize = 1;
    for cipherletter in key.iter() {
        if !cipherletter.is_empty() {
            ngram = cipherletter.len();
            break; 
        }
    }
    for (i, cipherletter) in key.iter().enumerate() {
        key_map.insert(
            cipherletter,
            &ALPHABET[i]
        );
    }
    for i in (0..ciphertext.chars().count()).step_by(ngram) {
        if key_map.contains_key(&ciphertext[i..i+ngram].to_string()) {
            plaintext += key_map[&ciphertext[i..i+ngram].to_string()];
        }
        else {
            plaintext += &"_".repeat(ngram);
        }
    }
    return String::from(plaintext);
}
    
fn score(plaintext: &str, ngram: &usize, ngram_log_p_table: &HashMap<String, f32>) -> f32 {
    let mut score: f32 = 0.0;
    for i in 0..(plaintext.chars().count() - ngram) {
        let ngram = &plaintext[i..i+ngram];
        if ngram_log_p_table.contains_key(ngram) {
            score += ngram_log_p_table[ngram];
        }
        else {
            score -= 8.0;
        }
    }
    return score;
}

fn main() {
    // Program Setup
    let args = Cli::from_args();
    let mut rng = thread_rng();

    // Clean english text sample source 
    let filter_alphanum_re = Regex::new(r"[^A-z1-9 ]").unwrap();
    let english_text_sample = include_str!("war_and_peace.txt");
    let english_text_stripped = filter_alphanum_re.replace_all(english_text_sample, "").to_uppercase();

    // Generate n-gram frequency table and compute logarithmic probability table for ngrams
    let mut ngram_freq_table = HashMap::<String, u32>::new();
    let mut total_ngrams = 0;
    for i in 0..(english_text_stripped.chars().count() - args.n_gram) {
        let ngram = &english_text_stripped[i..i+args.n_gram];
        if ngram.contains(" ") {
            continue;
        }
        if ngram_freq_table.contains_key(&ngram.to_string()) {
            ngram_freq_table.insert(
                ngram.to_string(),
                ngram_freq_table[&ngram.to_string()] + 1
            ); 
        }
        else {
            ngram_freq_table.insert(
                ngram.to_string(),
                1
            );
        }
        total_ngrams += 1;
    }
    let mut ngram_log_p_table = HashMap::<String, f32>::new();
    for key in ngram_freq_table.keys() {
        ngram_log_p_table.insert(
            key.to_string(),
            (ngram_freq_table[key] as f32 / total_ngrams as f32).log10()
        );
    }
    loop {
        // Create initial hypothesis key from ciphertext alphabet and ssc-ngram parameter
        let mut hypothesis_key: [String;26] = Default::default();
        let mut crypto_ngrams = HashSet::new();
        let cryptogram_u = args.cryptogram.to_uppercase();
        for i in (0..cryptogram_u.chars().count()).step_by(args.ngram_ssc) {
            crypto_ngrams.insert(&cryptogram_u[i..i+args.ngram_ssc]);
        }
        for (i, item) in crypto_ngrams.iter().enumerate() {
            hypothesis_key[i] = item.to_string();
        }
    // Search through key-space for most fit key according to ngram statistics
        let mut iterations: u32 = 0;
        let mut plaintext = build_plaintext(&args.cryptogram, &hypothesis_key);
        let mut fitness = score(&plaintext, &args.n_gram, &ngram_log_p_table);
        while iterations <= 1000 {
            let mut candidate_key = hypothesis_key.clone();
            candidate_key.swap(rng.gen_range(0,26),rng.gen_range(0,26));
            let candidate_plaintext = build_plaintext(&args.cryptogram, &candidate_key);
            let new_fitness = score(&candidate_plaintext, &args.n_gram, &ngram_log_p_table);
            if new_fitness > fitness {
                hypothesis_key = candidate_key.clone();
                fitness = new_fitness;
                plaintext = candidate_plaintext;
                iterations = 0;
            }
            else {
                iterations += 1;
            }
        }
        println!("\n{}\n{:?}", plaintext, hypothesis_key);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_plaintext_mono() {
        let ciphertext = String::from("ZYXWVUTSRQPONMLKJIHGFEDCBA");
        let key: [String;26] = [String::from("Z"),String::from("Y"),String::from("X"),
                                String::from("W"),String::from("V"),String::from("U"),
                                String::from("T"),String::from("S"),String::from("R"),
                                String::from("Q"),String::from("P"),String::from("O"),
                                String::from("N"),String::from("M"),String::from("L"),
                                String::from("K"),String::from("J"),String::from("I"),
                                String::from("H"),String::from("G"),String::from("F"),
                                String::from("E"),String::from("D"),String::from("C"),
                                String::from("B"),String::from("A")];
        assert_eq!(build_plaintext(&ciphertext, &key),
                   "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    }
    
    #[test]
    fn test_build_plaintext_digram() {
        let ciphertext = String::from("AABBCCDDEEFFGGHHIIJJKKLLMMNNOOPPQQRRSSTTUUVVWWXXYYZZ");
        let key: [String;26] = [String::from("AA"),String::from("BB"),String::from("CC"),
                                String::from("DD"),String::from("EE"),String::from("FF"),
                                String::from("GG"),String::from("HH"),String::from("II"),
                                String::from("JJ"),String::from("KK"),String::from("LL"),
                                String::from("MM"),String::from("NN"),String::from("OO"),
                                String::from("PP"),String::from("QQ"),String::from("RR"),
                                String::from("SS"),String::from("TT"),String::from("UU"),
                                String::from("VV"),String::from("WW"),String::from("XX"),
                                String::from("YY"),String::from("ZZ")];
        assert_eq!(build_plaintext(&ciphertext, &key),
                   "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    }
}
