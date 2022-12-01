// See https://aka.ms/new-console-template for more information
using rna_seq;

Console.WriteLine("Hello, World!");


var c = new Class1("./mRNA.csv");

c.Ttest();

c.PrintHistogram();

c.Filter();