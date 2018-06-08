library(coRdon)
context("codonTable class")

m <- matrix(sample(10:20, 128, replace = TRUE), ncol = 64)
df <- as.data.frame(m)
string1 <- "ATGGATTTTGGATTCTGTGACCTTCACAAACAGGCGTTGCCGGGCAAAAAGGCTTTG" # div. by 3
string2 <- "ATGCATGCAGTTGACCAGCTGTGACCTTCACACCGGGCAAAAAGGCTTGCATTGATAACG" # div. by 3
string3 <- "ATGGCGGGTATGACCATGCAGTTGACCAGCACCCACCATAATTTGCCTCATGGTAT" # not div. by 3
dna1 <- Biostrings::DNAStringSet(string1)
dna2 <- Biostrings::DNAStringSet(string2)
dna3 <- Biostrings::DNAStringSet(string3)
dna123 <- Biostrings::DNAStringSet(c(string1,string2, string3))
dna <- Biostrings::DNAStringSet(c(string1,string2))

test_that("codonTable works", {
    expect_message(codonTable(m), "alphabetically sorted codons")
    expect_equal(length(codonTable(m)), nrow(m))
    expect_equal(length(codonTable(rbind(m[1,]))), 1)
    expect_equal(length(codonTable(df)), nrow(df))
    expect_warning(codonTable(dna3), "not divisible by 3")
    expect_warning(codonTable(dna123), "not divisible by 3")
    expect_equal(getlen(codonTable(dna)), Biostrings::width(dna) %/% 3)
})

dna <- Biostrings::DNAStringSet(c(string1,string2, string1))
ct <- codonTable(dna)
ct <- setKO(ct, c("K00001", "K00002", "K00001"))
ct <- setCOG(ct, c("COG00001", "COG00002", "COG00001"))

test_that("subsetting works", {
    expect_equal(subset(ct, c(TRUE,FALSE,TRUE)), subset(ct, "K00001"))
    expect_equal(subset(ct, c(TRUE,FALSE,TRUE)), subset(ct, "COG00001"))
})
