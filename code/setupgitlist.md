Server Set-up
=============


## ssh



## git


## gitList


### apache & python


### customize the look of *gitList*
**Customize the content of the header bar.**



**Customize the style of text and other elements.** 
*gitList* is using TwitterBootstrap to make an automated static representation of the repositories. This is making use of `less`, an automated generator for cascading style sheets (css, tells the browser how to visualise certain elements of the website). The less files are stored in `var/www/gitlist/web/less/`. The created css file is stored in `var/www/gitlist/web/css`. So, The procedure for customization of style is 
 - find the `.less`-file that contains the element of question. For default text styles (headers, standard text, links) this should be `type.less`. For styles of inline code and code blocks this should be `code.less`. 
 - Modify the files. You need admin rights to do this, since the whole folder `www` is owned by the apache user `www-data`. So from the terminal call `sudo gedit filename.less`. 
 - update the `.css` file. in `var/www/gitlist/web/` you can just run `sudo make`, which executes the Makefile containing the proper prompt.  

Voila. That's it. 

#### less
As prerequisites for the customization, less needs to be properly installed in version 1.3.0 


##### node & npm

https://github.com/joyent/node/wiki/Installing-Node.js-via-package-manager#ubuntu

##### install less 
Don't install the recent version of less. Since some major changes made it incompatible with the mixins.less file in gitlist.

``` 
npm install -g less@1.3.0
```

