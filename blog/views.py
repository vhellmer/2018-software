from django.shortcuts import render

# Create your views here.
def front_page(request):
    return render(request, 'website/front_page.html', {})

def software(request):
	return render(request, 'website/software.html', {})

def discussion(request):
	return render(request, 'website/discussion.html', {})

def about(request):
	return render(request, 'website/about.html', {})