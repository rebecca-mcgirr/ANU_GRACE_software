      LOGICAL FUNCTION KBIT(IARRAY,IBIT)
C
C     KBIT
C
C 1.  KBIT PROGRAM SPECIFICATION
C
C 1.1.   KBIT is true if the IBIT-th bit in IARRAY is set (1),
C        false if it is not (0).
C        KBIT is designed to complement SBIT which sets or resets
C        bits identified in the same way.
C        Bits are numbered starting with 1 as the lowest-order bit
C        in the first word of IARRAY, 16 is the sign bit in the first
C        word of IARRAY, 17 is the low-order bit in the second word
C        of IARRAY, etc.
C
C 1.2.   RESTRICTIONS - NONE
C
C 1.3.   REFERENCES - FILE "SOLV2 = SOLV2 PURPOSE AND OVERALL STRUCTURE 
C 
C 2.  KBIT INTERFACE
C 
C 2.1.   CALLING SEQUENCE: CALL KBIT(IARRAY,IBIT) 
C 
C     INPUT VARIABLES:
C 
      integer*4 IARRAY(*) , ibit
C     - may or may not be an array in calling program.
C     - Variable in which the flag bits are located.
C 
C     IBIT = index of bit to test. Bits are numbered starting with
C            1 as the lowest order bit in IARRAY(1), 16 is the sign 
C            bit in IARRAY(1), 17 is the lowest order bit in IARRAY(2), 
C            etc. 
C 
C 
C     OUTPUT VARIABLES: 
C 
C     KBIT = FUNCTION VALUE = TRUE if the indicated bit is 1
C                           = FALSE if the indicated bit is 0 
C 
C 
C     CALLING SEQUENCE VARIABLE SPECIFICATION STATEMENTS: 
C 
C 
C 2.2.   COMMON BLOCKS USED 
C 
C 
      integer*4 MASK(32) 
C 
C     COMMON MSKCM CONTAINS THE STANDARD BIT MASK ARRAY 
C 
C     TO
C     VARIABLES IN MSKCM TO WHICH KBIT GIVES A VALUE: NONE
C 
C 
C     FROM
C     VARIABLES IN MSKCM FROM WHICH KBIT GETS A VALUE: MASK 
C     MASK(1)=1B, MASK(2)=2B ... MASK(16)=100000B (sign bit)
C     I.E. MASK(I) HAS I-TH BIT (SOLVE CONVENTION) SET
C 
C 
C 
C 
C 2.3.   DATA BASE ACCESSES 
C 
C 
C 2.4.   EXTERNAL INPUT/OUTPUT
C 
C 2.5.   SUBROUTINE INTERFACE:
C 
C     CALLING SUBROUTINES: THIS IS A SOLV2 UTILITY
C 
C     CALLED SUBROUTINES: 
C       UTILITIES IN SYSTEM: IAND(FUNCTION) 
C 
C 3.  LOCAL VARIABLES 
C 
C 
C END EQUIVALENCE 
C 
C 4.  CONSTANTS USED (description and references only. Data stmts below.) 
C 
C 
C 5.  INITIALIZED VARIABLES 
C 
C 
C 
C 
C 6.    LOCAL VARIABLES - SPECIFICATIONS AND DATA STATEMENTS

      integer*4 ia, ib
C 
C 6.1.  LOCAL VARIABLES - SPECIFICATION STATEMENTS
C 
C 6.2.  TEXT VARIABLES - SPECIFICATION,DATA STATEMENTS
C 
C 6.3.  INITIALIZED VARIABLES - DATA STATEMENTS 

      data mask /    1,      2,       4,       8,      16,    32,
     .              64,    128,     256,     512,    1024,  2048,
     .            4096,   8192,   16384,   32768,   65536,131072,
     .          262144,     524288,   1048576,  2097152, 4194304,
     .         8388608,   16777216,  33554432, 67108864,
     .       134217728,  268435456, 536870912, 
     .      1073741824,-2147483648 /
C 
C 
C 7.  PROGRAMMER: J.C.PIGG
C     LAST MODIFIED:
C# LAST COMPC'ED  820423:19:54 #
CMOD: JCP 1980 FEB 7: MOVED MASK FROM LOCAL TO COMMON MSKCM.
CMOD: JCP 79 JUL 27: CREATED
C 
C     PROGRAM STRUCTURE 
C 
C     1. Decompose IBIT into an array index IA and a bit index IB.
C 
      IA = (IBIT+31)/32 
      IB = IBIT - (IA-1)*32 
C
C
C     2. Test the appropriate bit
C
      KBIT = IAND(IARRAY(IA),MASK(IB)) .NE. 0
C
C
C     3. Return to calling program.
C
      END

