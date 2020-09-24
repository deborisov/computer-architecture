;-------------------Author: Borisov Dmitry cbse192------------
;-------------------Variant: 4--------------------------------
;-------------------Assigment: B[i] = -1 if A[i]<0------------
;-------------------           B[i] = 1 if A[i]>0-------------
;-------------------           B[i] = 0 if A[i]=0-------------
format PE console

entry start

include 'win32a.inc'

section '.data' data readable writable

        hHeap dd ?
        n dd ?
        A dd ?
        B dd ?
        tmpSize dd ?
        arrayPtr dd ?
        stackTop dd ?
        tmpElement dd ?

        inputNstr db 'Input positive size of array less than %d: ', 0
        notPositiveStr db 'Array size should be greater than 0.', 10, 13, 0
        greaterNumberStr db 'Array size should be less than %d.', 10, 13, 0
        outputStr db 'The output array is: ', 0
        endStr db 'Press any key to exit... ', 0
        inputElStr db 'Input %d-th element: ', 0
        endl db 10, 13, 0
        numberStr db '%d', 0
        spaceNumberStr db '%d ', 0


        NULL = 0
        MAX_ARRAY_SIZE = 1000; Can easily chnge this constant

section '.code' code readable executable

;---------------------------------Entry point of a programm-----------------------------------------------------------------
        start:
              call enterN
              call memAlloc
              call inputArray
              call formSignArray
              push [n]
              push [B]
              call printArray

              call freeHeap
              call endProgrammText
              push NULL
              call [ExitProcess]

;---------------------------------Procedure allocates memory in the heap for an arrays----------------------------------------
;---------------------------------Uses win32 API------------------------------------------------------------------------------
        memAlloc:
              invoke HeapCreate,0,0,0
              mov [hHeap], eax
              mov eax, [n]
              imul eax, 4
              invoke HeapAlloc, [hHeap], 0, eax
              mov [A], eax
              mov eax, [n]
              imul eax, 4
              invoke HeapAlloc, [hHeap], 0, eax
              mov [B], eax
              ret

;---------------------------------Procedure reads size of an array from the console-------------------------------------------
;---------------------------------And validate it, so it should be in the possible bounds-------------------------------------
        enterN:
              mov ebp, esp
              do:
                 cinvoke printf, inputNstr, MAX_ARRAY_SIZE
                 cinvoke scanf, numberStr, n
                 cmp [n], 0
                 jle less
                 cmp [n], MAX_ARRAY_SIZE
                 jae greater
                 jmp endWhile
              less:
                 cinvoke printf, notPositiveStr
                 jmp do
              greater:
                 cinvoke printf, greaterNumberStr, MAX_ARRAY_SIZE
                 jmp do
              endWhile:
                 mov esp, ebp
                 ret

;---------------------------------Procedure prints array-------------------------------------------
;---------------------------------Input: pointer to the start of an array, size of array-----------
;---------------------------------Input parameters should be pushed in the stack-------------------
        printArray:
                mov [stackTop], esp
                add esp, 4
                mov eax, [esp]
                mov [arrayPtr], eax
                add esp, 4
                mov eax, [esp]
                mov [tmpSize], eax
                mov esp, [stackTop]
                cinvoke printf, outputStr
                xor eax, eax
                printNumber:
                    mov ebx, eax
                    imul eax, 4
                    add eax, [arrayPtr]
                    cinvoke printf, spaceNumberStr, [eax]
                    mov eax, ebx
                    inc eax
                    cmp eax, [tmpSize]
                    jb printNumber
                cinvoke printf, endl
                ret

;-------------------------------------------Precedure calls getch and prints end programm text-----------------------------------
        endProgrammText:
                mov ebp, esp
                cinvoke printf, endStr
                call [getch]
                mov esp, ebp
                ret

;-------------------------------------------Free allocated memory-----------------------------------
;-------------------------------------------Hopefully that works------------------------------------
        freeHeap:
                invoke HeapFree, [hHeap], 0, [A]
                invoke HeapFree, [hHeap], 0, [B]
                invoke HeapDestroy, [hHeap]
                ret

;-------------------------------------------Procedure reads array of size n from the console-----------------------------------
        inputArray:
                mov ebp, esp
                xor eax, eax
                inputElement:
                        inc eax
                        mov ebx, eax
                        cinvoke printf, inputElStr, eax
                        mov eax, [A]
                        mov ecx, ebx
                        dec ecx
                        imul ecx, 4
                        add eax, ecx
                        cinvoke scanf, numberStr, eax
                        mov eax, ebx
                        cmp eax, [n]
                        jb inputElement
                mov esp, ebp
                ret

;-------------------------------------------Forms B array acording to the assigment specification-----------------------------------
        formSignArray:
                mov [stackTop], esp
                xor eax, eax
                formNextElement:
                        mov ebx, eax
                        mov ecx, eax
                        imul ebx, 4
                        imul ecx, 4
                        add ebx, [A]
                        add ecx, [B]
                        mov edx, [ebx]
                        cmp edx, 0
                        je zero
                        jl minus
                        jg plus
                        zero:
                                mov edx, 0
                                mov [ecx], edx
                                jmp   checkEndCycle
                        minus:
                                mov edx, -1
                                mov [ecx], edx
                                jmp checkEndCycle
                        plus:
                                mov edx, 1
                                mov [ecx], edx
                                jmp  checkEndCycle
                        checkEndCycle:
                                inc eax
                                cmp eax, [n]
                                jb formNextElement
                mov esp, [stackTop]
                ret

section '.idata' import data readable

        library kernel, 'kernel32.dll', \
                msvcrt, 'msvcrt.dll'

        import kernel, \
               ExitProcess, 'ExitProcess', \
               HeapCreate, 'HeapCreate',\
               HeapAlloc, 'HeapAlloc', \
               HeapFree, 'HeapFree', \
               HeapDestroy, 'HeapDestroy'

        import msvcrt,\
              printf, 'printf',\
              scanf, 'scanf',\
              getch, '_getch'