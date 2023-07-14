# PROGRAM:     time_extract.awk
# AUTHOR:      Anthony Purcell
# DATE:        October 3 2012
# DESCRIPTION: Takes a Date entry from the gamit doy program and extracts
#              hour and minute values.
#              When called the value of the variable 'string' should be set on 
#              the command line.
#              For string == 0 the values are output as numerical values.
#              For string == 1 the values are output as two character strings.

{
 num = split($3, time, ":")
 if ( string == 0 ) { printf("%g %g\n", time[1], time[2]) }
 if ( string == 1 )
   {
    if ( length(time[1]) == 1 ) { time[1] = "0" time[1] }
    printf("%-s %-s\n", time[1], time[2])
   }
}
