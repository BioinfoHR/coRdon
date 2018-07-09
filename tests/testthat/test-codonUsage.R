library(coRdon)
context("codonUsage")

string1 <- "ATGGATTTTGGATTCTGTGACCTTCACAAACAGGCGTTGCCGGGCAAAAAGGCTTTG"
string2 <- "ATGCATGCAGTTGACCAGCTGTGACCTTCACACCGGGCAAAAAGGCTTGCATTGATAACG"
dna <- Biostrings::DNAStringSet(c(string1,string2))
ct <- codonTable(dna)

test_that("genCode works", {
    expect_equal(dim(MILC(ct)), dim(MILC(ct, id_or_name2 = "2")))
})

test_that("CU methods produce properly structured output", {
    expect_equal(dim(MILC(ct)), dim(B(ct)))
    expect_equal(dim(MILC(ct)), dim(MCB(ct)))
    expect_equal(dim(MILC(ct)), dim(ENCprime(ct)))
    expect_equal(dim(ENC(ct)), dim(SCUO(ct)))
})

test_that("expressivity statistics produce properly structured output", {
    expect_equal(dim(MELP(HD59, ribosomal = TRUE)),
                 dim(E(HD59, ribosomal = TRUE)))
    expect_equal(dim(MELP(HD59, ribosomal = TRUE)),
                 dim(CAI(HD59, ribosomal = TRUE)))
    expect_equal(dim(MELP(HD59, ribosomal = TRUE)),
                 dim(Fop(HD59, ribosomal = TRUE)))
})

test_that("CU calculations produce no NAs/NaNs", {
    expect_false(any(is.na(MILC(HD59, ribosomal = TRUE))))
    expect_false(any(is.na(MELP(HD59, ribosomal = TRUE))))
    expect_false(any(is.na(B(HD59, ribosomal = TRUE))))
    expect_false(any(is.na(E(HD59, ribosomal = TRUE))))
    expect_false(any(is.na(CAI(HD59, ribosomal = TRUE))))
})
