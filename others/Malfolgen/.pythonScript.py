#!/usr/bin/env python

import random
import time

def malfolge(z):
    l = range(1,11)
    random.shuffle(l)
    while l:
        #print l
        print '\n*** Nur noch ' + str(len(l)) + ' Aufgabe(n) in der ' + str(z) +'er Malfolge. ***\n'
        a = l[-1]
        s = str(a) + ' mal ' + str(z) + ' = '
        eingabe = int(raw_input(s))
        if eingabe == a * z:
            print '\nRichtig! :)'
            l.pop()
        else:
            print '\nLeider falsch. :( Merke dir: ' + s + str(a*z)
            l.append(a)
            random.shuffle(l)


while True:
    print '\n'
    print '[b]: Moechtest du eine [b]estimmte Malfolge lernen?'
    print '[z]: Moechtest du eine [z]ufaellig ausgewaehlte Malfolge lernen?'
    print '[m]: Moechtest du gleich [m]ehrere Malfolgen lernen?'
    print '[a]: Moechtest du [a]ufhoeren zu lernen?'
    print '\n'
    eingabe = raw_input()[0].lower()
    if eingabe == 'b':
        eingabe = int(raw_input('Welche Malfolge Moechtest du lernen? '))
        malfolge(eingabe)
    elif eingabe == 'z':
        eingabe = int(raw_input('Bis zur welcher Malfolge hast du schon gelernt? '))
        malfolge(random.randint(2, eingabe))
    elif eingabe == 'm':
        eingabe = int(raw_input('Bis zur welcher Malfolge hast du schon gelernt? '))
        l = range(2, eingabe+1)
        random.shuffle(l)
        while l:
            print '\n****** Nur noch ' + str(len(l)) + ' Malfolge(n) zu ueben. ******\n'
            malfolge(l.pop())
    elif eingabe == 'a':
        print 'Auf Wiedersehen, hat Spass gemacht mit dir. :)'
        time.sleep(4)
        break
    else:
        print 'Diese Eingabe kenne ich nicht. :('
