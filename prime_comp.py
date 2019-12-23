#!/usr/bin/python


def reverse_string(s):
    """Return a reversed copy of `s`"""
    return s[::-1]

#the script calculates the reverse compliment of any given sequence espicially primers
#merit:: if sequencing length is longer than your amplicon then you are sequencing the reversed (reverse primer in the end of the F strand)


For_Primer = "GTGCCAGCMGCCGCGGTAA"
Rev_Primer = "GGACTACHVGGGTWTCTAAT"


Rev_For = reverse_string(For_Primer)
print("the Rev_For_primer is " + Rev_For)

Rev_Rev = reverse_string(Rev_Primer)
print("the Rev_Rev_primer is " + Rev_Rev)

Rev_For_Primer_A_replaced = Rev_For.replace("A", "t")
Rev_For_Primer_T_replaced = Rev_For_Primer_A_replaced.replace("T", "a")
Rev_For_Primer_C_replaced = Rev_For_Primer_T_replaced.replace("C", "g")
Rev_For_Primer_G_replaced = Rev_For_Primer_C_replaced.replace("G", "c")
Rev_For_Primer_Y_replaced = Rev_For_Primer_G_replaced.replace("Y", "r")
Rev_For_Primer_R_replaced = Rev_For_Primer_Y_replaced.replace("R", "y")
Rev_For_Primer_K_replaced = Rev_For_Primer_R_replaced.replace("K", "m")
Rev_For_Primer_M_replaced = Rev_For_Primer_K_replaced.replace("M", "k")
Rev_For_Primer_B_replaced = Rev_For_Primer_M_replaced.replace("B", "v")
Rev_For_Primer_D_replaced = Rev_For_Primer_B_replaced.replace("D", "h")
Rev_For_Primer_H_replaced = Rev_For_Primer_D_replaced.replace("H", "d")
Rev_For_Primer_V_replaced = Rev_For_Primer_H_replaced.replace("V", "b")

Rev_For_Primer_comp = Rev_For_Primer_V_replaced

print("The reverse compliment sequence of the forward primer to be added to the end of Reverse primer is "+ Rev_For_Primer_comp.upper())

Rev_Rev_Primer_A_replaced = Rev_Rev.replace("A", "t")
Rev_Rev_Primer_T_replaced = Rev_Rev_Primer_A_replaced.replace("T", "a")
Rev_Rev_Primer_C_replaced = Rev_Rev_Primer_T_replaced.replace("C", "g")
Rev_Rev_Primer_G_replaced = Rev_Rev_Primer_C_replaced.replace("G", "c")
Rev_Rev_Primer_Y_replaced = Rev_Rev_Primer_G_replaced.replace("Y", "r")
Rev_Rev_Primer_R_replaced = Rev_Rev_Primer_Y_replaced.replace("R", "y")
Rev_Rev_Primer_K_replaced = Rev_Rev_Primer_R_replaced.replace("K", "m")
Rev_Rev_Primer_M_replaced = Rev_Rev_Primer_K_replaced.replace("M", "k")
Rev_Rev_Primer_B_replaced = Rev_Rev_Primer_M_replaced.replace("B", "v")
Rev_Rev_Primer_D_replaced = Rev_Rev_Primer_B_replaced.replace("D", "h")
Rev_Rev_Primer_H_replaced = Rev_Rev_Primer_D_replaced.replace("H", "d")
Rev_Rev_Primer_V_replaced = Rev_Rev_Primer_H_replaced.replace("V", "b")

Rev_Rev_Primer_comp = Rev_Rev_Primer_V_replaced

print("The reverse compliment sequence of the reverse primer to be added to the end of forward primer is "+ Rev_Rev_Primer_comp.upper())
