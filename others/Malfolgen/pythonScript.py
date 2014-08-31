#!/usr/bin/env python

import random
import time
from pylab import array

def malfolge(z):
    nAufgaben = 0
    nFehler = 0
    l = range(1,11)
    random.shuffle(l)
    while l:
        nAufgaben += 1
        #print l
        print '\n*** Nur noch ' + str(len(l)) + ' Aufgabe(n) in der ' + str(z) +'er Malfolge. ***\n'
        a = l[-1]
        s = str(a) + ' mal ' + str(z) + ' = '
        eingabe = int(raw_input(s))
        if eingabe == a * z:
            print '\nRichtig! :)'
            l.pop()
        elif eingabe == a * z * 100:
            print '\nRichtig! :) (Aufgabe bleibt aber im Test)'
            random.shuffle(l)
        else:
            nFehler += 1
            print '\nLeider falsch. :( Merke dir: ' + s + str(a*z)
            l.append(a)
            random.shuffle(l)
    return nAufgaben - nFehler , nAufgaben

def printSchnitt(schnitt, v=True):
    nRichtig = schnitt[0]
    nAufgaben = schnitt[1]
    prozent = float(nRichtig) / nAufgaben * 100.
    print '\n +++ Du hast ' + str(nRichtig) + ' von ' + str(nAufgaben) + ' Aufgaben Richtig. (' + str(prozent) + ' Prozent) +++ \n'

schnitt = array([0,0])
while True:
    print '\n'
    print '[b]: Moechtest du eine [b]estimmte Malfolge lernen?'
    print '[z]: Moechtest du eine [z]ufaellig ausgewaehlte Malfolge lernen?'
    print '[m]: Moechtest du gleich [m]ehrere Malfolgen lernen?'
    print '[a]: Moechtest du [a]ufhoeren zu lernen?'
    print '\n(Wenn du zwei Nullen an die Zahl haengst, dann bleibt die Aufgabe im Test. (z.B. 1500 statt 15))'
    print '\n'
    eingabe = raw_input()[0].lower()
    if eingabe == 'b':
        eingabe = int(raw_input('Welche Malfolge Moechtest du lernen? '))
        schnitt += malfolge(eingabe)
        printSchnitt(schnitt)
    elif eingabe == 'z':
        eingabe = int(raw_input('Bis zur welcher Malfolge hast du schon gelernt? '))
        schnitt += malfolge(random.randint(2, eingabe))
        printSchnitt(schnitt)
    elif eingabe == 'm':
        eingabe = int(raw_input('Bis zur welcher Malfolge hast du schon gelernt? '))
        l = range(2, eingabe+1)
        random.shuffle(l)
        while l:
            print '\n****** Nur noch ' + str(len(l)) + ' Malfolge(n) zu ueben. ******\n'
            schnitt += malfolge(l.pop())
            printSchnitt(schnitt)
    elif eingabe == 'a':
        print 'Auf Wiedersehen, hat Spass gemacht mit dir. :)'
        time.sleep(4)
        break
    else:
        print 'Diese Eingabe kenne ich nicht. :('
