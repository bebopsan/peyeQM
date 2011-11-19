# -*- coding: utf-8 -*-
## module inform
"""
Created on Fri Nov 18 21:10:37 2011

@author: santiago
"""
def inform(text):
    import tweepy
    
    # == OAuth Authentication ==
    #
    # This mode of authentication is the new preferred way
    # of authenticating with Twitter.
    
    # The consumer keys can be found on your application's Details
    # page located at https://dev.twitter.com/apps (under "OAuth settings")
    consumer_key="YBva1alKSOBtHZfRgOo41Q"
    consumer_secret="gERjZ97u3g9hj67fFFtnIDUqr2bzmu5HVxP0Tgzs"
    
    # The access tokens can be found on your applications's Details
    # page located at https://dev.twitter.com/apps (located
    # under "Your access token")
    access_token="46562492-SZW500iBuNch7TN9WsUuyehmf1PCOdGhpbuAM8fiJ"
    access_token_secret="uQlNY4IpUZcOuXG7ezRzIpCTwTBa9Yk64AYcCd3M"
    
    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
    auth.set_access_token(access_token, access_token_secret)
    
    api = tweepy.API(auth)
    
    # If the authentication was successful, you should
    # see the name of the account print out
    print api.me().name
    
    # If the application settings are set for "Read and Write" then
    # this line should tweet out the message to your account's
    # timeline. The "Read and Write" setting is on https://dev.twitter.com/apps
    #api.update_status('@Bebopsan OAuth authentication via Tweepy!')
    me = api.me()
    help(me)
    api.send_direct_message(user_id = me.id, text = text)


