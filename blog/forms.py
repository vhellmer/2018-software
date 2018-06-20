from django import forms

class NameForm(forms.Form):
    input1 = forms.IntegerField(label='input1')