;-------------------Author: Borisov Dmitry cbse192------------
;-------------------Variant: 4--------------------------------
;-------------------Assigment: count number of vovels-------------
;-------------------and consonants in the input string-------------
format PE console

entry start

include 'win32a.inc'

section '.data' data readable writable
        endStr db 'Press any key to exit... ', 0;prints in the end of program
        stringForInput db 'Input your string: ', 0;prints before user inputs string
        endl db 10, 13, 0;end of line
        numberStr db '%d', 0
        stringStr db '%s', 0
        inputStr db 1001 dup (?), 0;User input
        numberVovelsStr db 'Number of vowels: %d', 0;Message for number of vovels in string
        numbrerConsonantsStr db 'Number of consonants: %d', 0;Message for number of consonants in string
        numberVovels dd 0;mark for number of vowels
        numberConsonants dd 0;mark for number of consonants
        NULL = 0

section '.code' code readable executable

;---------------------------------Entry point of a programm-----------------------------------------------------------------
        start:
              call readInput

              stdcall countAnswer, inputStr;Main program logic
              add esp, 4
              call printAnswer

              call endProgrammText;Call getch and after close down program
              push NULL
              call [ExitProcess]

        ;------------------------------------------------------------------
        ; Input: string to process as stack argument
        ; Procedure counts number of vovels in numberVovels
        ; and number of consonants in numberConsonants
        countAnswer:
                mov esi, [esp + 4];using input string as stack argument
                cld;direction from begin to end
                xor ecx, ecx;End loop condition char equals to 0 (ecx)
                loopStr:
                        lodsb;load char
                        cmp cl, al;if zero goto end of loop
                        je endSrc
                        cmp al, 'A';Less than 'A' - not a letter
                        jl loopStr;Read next char
                        cmp al, 'Z';Less than 'Z' and greater than 'A' - greater letter
                        jle greaterChar;Precess greater letter
                        cmp al, 'z';Greater than 'z' - not a letter
                        jg loopStr;Read next char
                        cmp al, 'a';Less than 'z' and greater than 'a' - lower letter
                        jge lowerChar;Precess lower letter
                        jmp loopStr;Read next char
                lowerChar:
                        ;If any vovel increase number of vovels
                        cmp al, 'a'
                        je incVovel
                        cmp al, 'e'
                        je incVovel
                        cmp al, 'i'
                        je incVovel
                        cmp al, 'o'
                        je incVovel
                        cmp al, 'u'
                        je incVovel
                        cmp al, 'y'
                        je incVovel

                        ;If not increase number of consonants
                        jmp incConsonants
                greaterChar:
                        ;If any vovel increase number of vovels
                        cmp al, 'A'
                        je incVovel
                        cmp al, 'E'
                        je incVovel
                        cmp al, 'I'
                        je incVovel
                        cmp al, 'O'
                        je incVovel
                        cmp al, 'U'
                        je incVovel
                        cmp al, 'Y'
                        je incVovel

                        ;If not increase number of consonants
                        jmp incConsonants
                incVovel:
                        inc [numberVovels];Increase number of vovels
                        jmp loopStr;Read next char
                incConsonants:
                        inc [numberConsonants];Increase number of consonants
                        jmp loopStr;Read next char
                endSrc:
                        ret;Return from procedure

        ;------------------------------------------------------------------
        ; Procedure prints two strings with number of vovels in first
        ; and number of consonants in second
        printAnswer:
                cinvoke printf, numberVovelsStr, [numberVovels];print number of vovels
                cinvoke printf, endl;print new string
                cinvoke printf, numbrerConsonantsStr, [numberConsonants];print number of consonants
                cinvoke printf, endl;print new string
                ret

        ;------------------------------------------------------------------
        ; Procedure prints text before the program will
        ; close down and calls getch
        endProgrammText:
                cinvoke printf, endStr
                call [getch]
                ret

        ;------------------------------------------------------------------
        ; Procedure prints prompt for enter string and reads string into
        ; strignForInput mark
        readInput:
                cinvoke printf, stringForInput
                cinvoke scanf, stringStr, inputStr
                ret


section '.idata' import data readable

        library kernel, 'kernel32.dll', \
                msvcrt, 'msvcrt.dll'

        import kernel, \
               ExitProcess, 'ExitProcess'

        import msvcrt,\
              printf, 'printf',\
              scanf, 'scanf',\
              getch, '_getch'