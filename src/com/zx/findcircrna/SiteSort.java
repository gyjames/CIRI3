package com.zx.findcircrna;

public class SiteSort implements Comparable<SiteSort> {
	private Integer site;
	private String[] length;
   
 
    public SiteSort(int site,String[] length) {
    	this.site = site;
    	this.length = length;    
       }
   
    public String[] getLength() {
		return length;
	}




	public Integer getSite() {
		return site;
	}

	@Override
    public int compareTo(SiteSort siteSort) {
        return this.site-siteSort.getSite();     
    }
	
	}

