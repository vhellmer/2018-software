from django.shortcuts import render
from django.http import HttpResponseRedirect
from .forms import NameForm


# Create your views here.
def front_page(request):
    return render(request, 'website/front_page.html', {})

def software(request):
	return render(request, 'website/software.html', {})

def discussion(request):
	return render(request, 'website/discussion.html', {})

def about(request):
	return render(request, 'website/about.html', {})

def addition(request):
	try:
		sum1 = request.GET['input1']
	except ValueError:
		sum1 = -1
	return sum1

def inputs(request):
    # if this is a POST request we need to process the form data
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = NameForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            return HttpResponseRedirect('/software/')

    # if a GET (or any other method) we'll create a blank form
    else:
        form = NameForm()

    return render(request, 'website/front_page.html', {'form': form})