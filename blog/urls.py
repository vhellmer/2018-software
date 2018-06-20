from django.conf.urls import url
from django.urls import include, path

from . import views
urlpatterns = [
    url(r'^$', views.front_page, name='front_page'),
    url(r'^software/$', views.software, name='software'),
    url(r'^discussion/$', views.discussion, name='discussion'),
    url(r'^about/$', views.about, name='about'),
    path('do_stuff/', views.inputs, name='inputs'),
]