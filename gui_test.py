from Tkinter import *
import tkMessageBox

#Create window
root = Tk()
root.title("O-ribosome Software Tool")
root.geometry("500x500")
root.configure(background="slategray")


#Makes input frame
input_frame = Frame(root)
input_frame.pack()

#Makes submit frame
submit_frame = Frame(root)
submit_frame.pack(side = BOTTOM)

#Makes input box
L1 = Label(input_frame, text="Expression system", font=('arial', 16))
L1.pack( side = LEFT)
E1 = Entry(input_frame, bd =5, relief=FLAT, highlightbackground='darkgray', font=('arial', 16))
E1.pack(side = RIGHT)

#Message pop-up
def helloCallBack():
   tkMessageBox.showinfo( "Hello Python", E1.get())


#Makes button 
B = Button(submit_frame, text ="Submit", highlightbackground='slategray', command = helloCallBack, font=('arial', 16))
B.pack(side = LEFT)


root.mainloop()



